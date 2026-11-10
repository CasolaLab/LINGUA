#!/usr/bin/env python3
"""
run_monitor.py

Watches a long run on Dayhoff and the machine's health, pausing the run (SIGSTOP to its process group, nothing is
lost) when something looks unsafe and resuming it (SIGCONT) when it is safe again. Emails on every pause, resume,
unreadable sensor and at the end of the run. Logs one line per interval to <log>.

Pause when any of:
  CPU Tctl >= --pause-c for --pause-n consecutive samples
  any fan reads 0 rpm
  free space on the run's filesystem < --min-free-gb
  the NVMe link watcher's last line shows a link other than the baseline, or kernel error lines > 0
Resume when all are clear and Tctl <= --resume-c for --resume-n consecutive samples.
Sensors that cannot be read are reported ("unknown, look at this"), never acted on.
Alert emails (pause, resume, unreadable sensor) are sent at most once per --mail-min-interval seconds; suppressed
alerts are counted in the next one. The start and end emails are always sent. At start, processes of the group found
stopped (state T) are resumed with SIGCONT, so restarting the monitor can never leave the run paused.
The run is finished when <run_root>/STATUS no longer says RUNNING; the monitor then emails that status and the
validation totals and exits. It never kills the run.

Usage: python3 run_monitor.py <run_root> <pgid> <log> [--interval 60] [--pause-c 88] [--resume-c 78]
                              [--mail-min-interval 1800]
"""
import argparse
import datetime
import glob
import os
import shutil
import signal
import smtplib
import time
from email.message import EmailMessage

RECIPIENT = "ccasola@tamu.edu"


def mail(subject, body):
    try:
        m = EmailMessage()
        m["From"] = m["To"] = RECIPIENT
        m["Subject"] = "Dayhoff: " + subject
        m.set_content(body)
        with smtplib.SMTP("smtp-relay.tamu.edu", 25, timeout=60) as s:
            s.send_message(m)
        return True
    except Exception as e:
        return f"mail failed: {e}"


def hwmon(name):
    for h in glob.glob("/sys/class/hwmon/hwmon*"):
        try:
            if open(h + "/name").read().strip() == name:
                return h
        except OSError:
            pass
    return None


def read_int(p):
    try:
        return int(open(p).read().strip())
    except (OSError, ValueError):
        return None


def link_watch_last():
    files = sorted(glob.glob(os.path.expanduser("~/link_watch_*.csv")), key=os.path.getmtime)
    if not files:
        return None
    lines = open(files[-1]).read().strip().split("\n")
    return dict(zip(lines[0].split(","), lines[-1].split(","))) if len(lines) > 1 else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_root")
    ap.add_argument("pgid", type=int)
    ap.add_argument("log")
    ap.add_argument("--interval", type=int, default=60)
    ap.add_argument("--pause-c", type=float, default=88)
    ap.add_argument("--pause-n", type=int, default=3)
    ap.add_argument("--resume-c", type=float, default=78)
    ap.add_argument("--resume-n", type=int, default=5)
    ap.add_argument("--min-free-gb", type=float, default=100)
    ap.add_argument("--mail-min-interval", type=int, default=1800)
    a = ap.parse_args()
    last_alert, suppressed = 0.0, []

    def alert(subject, body):
        nonlocal last_alert, suppressed
        if time.time() - last_alert < a.mail_min_interval:
            suppressed.append(f"{datetime.datetime.now():%H:%M} {subject}")
            return
        if suppressed:
            body += "\n\nSuppressed since the last alert:\n" + "\n".join(suppressed)
            suppressed = []
        mail(subject, body)
        last_alert = time.time()

    stopped = []
    for d in glob.glob("/proc/[0-9]*"):
        try:
            st = open(d + "/stat").read().rsplit(")", 1)[1].split()
            if int(st[2]) == a.pgid and st[0] == "T":
                stopped.append(d)
        except (OSError, IndexError, ValueError):
            pass
    if stopped:
        os.killpg(a.pgid, signal.SIGCONT)
    k10, smm = hwmon("k10temp"), hwmon("dell_smm")
    base = link_watch_last()
    base_link = base.get("ssd0_link") if base else None
    new = not os.path.exists(a.log)
    log = open(a.log, "a", buffering=1)
    if new:
        log.write("time,tctl_c,max_tccd_c,fan1,fan2,fan3,load1,free_gb,ssd0_link,kernel_err_lines,run_status,paused,note\n")
    paused, hot, cool, warned_unknown = False, 0, 0, False
    mail("stage 3 full run monitor started", f"Run: {a.run_root}\nProcess group: {a.pgid}\n"
         f"Pause at Tctl >= {a.pause_c} C x{a.pause_n}, fan 0 rpm, free < {a.min_free_gb} GB, NVMe link change or "
         f"kernel errors; resume at Tctl <= {a.resume_c} C x{a.resume_n}.\nAlerts at most every {a.mail_min_interval // 60} min.\n"
         f"Stopped processes found and resumed at start: {len(stopped)}\nBaseline NVMe link: {base_link}\nLog: {a.log}")
    while True:
        tctl = read_int(k10 + "/temp1_input") if k10 else None
        tccd = [read_int(p) for p in glob.glob(k10 + "/temp*_input")] if k10 else []
        tccd = [x for x in tccd if x is not None]
        fans = [read_int(f"{smm}/fan{i}_input") if smm else None for i in (1, 2, 3)]
        load1 = os.getloadavg()[0]
        free = shutil.disk_usage(a.run_root).free / 1e9
        lw = link_watch_last()
        link = lw.get("ssd0_link") if lw else None
        kerr = lw.get("kernel_err_lines") if lw else None
        try:
            status = open(os.path.join(a.run_root, "STATUS")).read().strip()
        except OSError:
            status = "NO STATUS"
        tc = tctl / 1000 if tctl is not None else None
        reasons, unknown = [], []
        if tc is None:
            unknown.append("CPU temperature")
        else:
            hot = hot + 1 if tc >= a.pause_c else 0
            cool = cool + 1 if tc <= a.resume_c else 0
            if hot >= a.pause_n:
                reasons.append(f"Tctl {tc:.1f} C for {hot} samples")
        if any(f is None for f in fans):
            unknown.append("fan speed")
        elif any(f == 0 for f in fans):
            reasons.append(f"fan at 0 rpm {fans}")
        if free < a.min_free_gb:
            reasons.append(f"only {free:.0f} GB free")
        if lw is None:
            unknown.append("NVMe link watcher log")
        else:
            if base_link and link != base_link:
                reasons.append(f"NVMe link now '{link}' (was '{base_link}')")
            if kerr not in (None, "0"):
                reasons.append(f"kernel error lines {kerr}")
        note = ""
        running = status.startswith("RUNNING")
        if running and reasons and not paused:
            os.killpg(a.pgid, signal.SIGSTOP)
            paused, note = True, "PAUSED: " + "; ".join(reasons)
            alert("stage 3 run PAUSED", note + f"\nRun: {a.run_root}\nResume by hand: kill -CONT -{a.pgid}")
        elif paused and not reasons and cool >= a.resume_n:
            os.killpg(a.pgid, signal.SIGCONT)
            paused, note = False, "RESUMED"
            alert("stage 3 run RESUMED", f"Conditions clear (Tctl {tc:.1f} C).\nRun: {a.run_root}")
        if unknown and not warned_unknown:
            warned_unknown = True
            alert("stage 3 monitor: sensor unreadable, look at this", "Cannot read: " + ", ".join(unknown) +
                  "\nNo action taken.")
        log.write(",".join(str(x) for x in (
            datetime.datetime.now().isoformat(timespec="seconds"), f"{tc:.1f}" if tc is not None else "NA",
            f"{max(tccd) / 1000:.1f}" if tccd else "NA", *[f if f is not None else "NA" for f in fans],
            f"{load1:.2f}", f"{free:.0f}", (link or "NA").replace(",", ";"), kerr or "NA", status.replace(",", ";"),
            int(paused), note.replace(",", ";"))) + "\n")
        if not running and status != "NO STATUS":
            body = f"Status: {status}\nRun: {a.run_root}\n"
            vs = os.path.join(a.run_root, "validation_summary.tsv")
            if os.path.exists(vs):
                rows = [l.rstrip("\n").split("\t") for l in open(vs)][1:]
                t = [r for r in rows if r[1] != "alignments"]
                body += f"Tables byte-identical: {sum(r[3] == 'BYTE-IDENTICAL' for r in t)} of {len(t)}\n"
                body += "\n".join(f"{r[0]}: {r[3]}" for r in rows if r[1] == "alignments") + "\n"
            if suppressed:
                body += "Suppressed alerts:\n" + "\n".join(suppressed) + "\n"
            body += f"Monitor log: {a.log}"
            mail(f"stage 3 full run {status.split()[0]}", body)
            break
        time.sleep(a.interval)


if __name__ == "__main__":
    main()
