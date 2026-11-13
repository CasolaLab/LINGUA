"""
make_workflows.py

Workflow figures of the six LINGUA stages (0-5), without counts from any dataset, one PNG and one SVG per stage,
plus one PowerPoint deck with one stage per slide. The step texts follow each stage README; thresholds are the
defaults of the code. Figures are drawn at slide size so every text is at least 14 pt on the slide.

Usage (from the LINGUA root):
    python docs/figures/make_workflows.py [--overwrite]
Output: docs/figures/stage<N>_workflow.{png,svg} and docs/figures/LINGUA_workflows.pptx
Requires matplotlib and python-pptx.
"""
import os
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUTDIR = os.path.join(os.getcwd(), "docs", "figures")
OVERWRITE = "--overwrite" in sys.argv
mpl.rcParams.update({"font.family": "Arial", "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.cal": "Arial", "mathtext.it": "Arial", "mathtext.bf": "Arial:bold", "svg.fonttype": "none", "pdf.fonttype": 42})
HUE, INK, NOTE, IOFILL = "#14694c", "#1a1a1a", "#5b5b5b", "#eef3f1"
FS_T, FS_B = 15, 14
MIN_FS = 14
FIG_W, MAX_H = 12.9, 6.5
LINE = 0.245
PAD = 0.12
GAP_X, GAP_Y = 0.26, 0.55


def io(title, body):
    return dict(kind="io", title=title, body=body, note=None)


def step(title, body, note=None, optional=False):
    return dict(kind="opt" if optional else "step", title=title, body=body, note=note)


STAGES = [
    (0, "Protein preprocessing", [[
        io("Input", "GFF3 and protein FASTA\nper species, chosen\nupstream (for example\nby BUSCO)"),
        step("Protein QC", "removes proteins with\nno start methionine,\nshorter than 30 aa, or\nwith an internal stop"),
        step("Isoform selection", "one protein per gene\n(the longest isoform);\ngene IDs standardized"),
        step("CDS and cDNA", "longest CDS and cDNA\nof each kept isoform", note="optional", optional=True),
        io("Output", "one standardized,\nnon-redundant proteome\nper species\n→ stage 1"),
    ]], None),
    (1, "Gene stratigraphy", [[
        io("Input", "standardized proteomes\n(stage 0); rooted species\ntree (optional)"),
        step("OrthoFinder", "DIAMOND all against all;\nhierarchical orthogroups\n(HOGs) and a node-labeled\nspecies tree"),
        step("Node placement", "each gene placed at the\nnode where its HOG first\nappears; unclustered genes\nare species-specific",
             note="HOG table of the clade root"),
    ], [
        step("Cross-check", "removes candidates whose\northogroup or pairwise\northologues include\nanother species"),
        step("Verification", "re-derives every candidate\nfrom the HOG table; stops\non any disagreement", note="check, removes no gene"),
        io("Output", "candidate lineage-specific\ngenes (cLSGs), protein\nFASTA per species\n→ stage 2"),
    ]], None),
    (2, "Homolog exclusion", [
        [io("Input", "candidate genes (stage 1), one protein FASTA per species")],
        [step("InterProScan", "member databases;\nPANTHER family\nrule (optional)"),
         step("CDD", "RPS-BLAST\nE ≤ 0.01, ≥ 50%\nof the domain"),
         step("BLASTp", "database outside\nthe clade; E ≤ 10$^{-3}$,\n≥ 50% of query"),
         step("jackhmmer", "3 iterations;\nE ≤ 10$^{-3}$, ≥ 50%\nof query"),
         step("TE proteins", "BLASTp vs a TE\nprotein library,\n≥ 50% of the gene"),
         step("Organelle", "tBLASTn vs\norganelle genomes;\nany hit counts")],
        [step("Combine", "a qualifying hit in any\nscreen removes the gene;\nan unfinished search\nnever passes"),
         step("HOG rule", "removes genes sharing\n≥ 20 aa with the hit region\nof a removed gene\nof their HOG"),
         step("Paralog rule", "removes genes with ≥ 50%\nidentity over ≥ 50% of\nlength to a non-candidate\ngene of their proteome"),
         io("Output", "lineage-specific genes,\nprotein FASTA and BED of\nCDS coordinates\n→ stage 3")],
    ], "Screens run on all genes or each on the survivors of the previous one. Hits to the analysis species are not counted.\n"
       "Thresholds shown are the defaults; each can be changed."),
    (3, "Alignment extraction", [[
        io("Input", "HAL alignment\n(Progressive Cactus);\nreference BED and CDS\nFASTA; species tree"),
        step("1  Extract blocks", "mafExtractor at each\nexon; stops on a failed\nor empty extraction"),
        step("2–3  Dominant chain", "each species keeps its\ncollinear block chain\nper exon; jumps > 50 kb,\nchromosome or strand\nare dropped"),
        step("4–7  Exon FASTA", "blocks joined into an\naligned FASTA, strand\ncorrected and validated"),
        step("8–8c  Stitch exons", "one row per species and\nancestor from its\ncollinear exon chain\n(≤ 50 kb); missing-data\nmask written"),
    ], [
        step("9–11  Trim and check", "extra base per exon\nremoved; the reference\nmust equal the CDS"),
        step("12–14b  Mutations", "start and stop codons;\nframeshift deletions\n(masked gaps skipped);\nevidence per deletion"),
        step("15–16  Motifs", "first and last three\nbases; start and stop\nmotifs checked"),
        step("17  Combine", "per-species span of\nthe kept exons, extant\nor ancestral"),
        io("Output", "alignment per locus;\nindel and stop tables;\ncoordinates\n→ stages 4 and 5"),
    ]], None),
    (4, "Enablers", [[
        io("Input", "stage 3 results;\ncandidate proteins\n(stage 2);\ntree_config.json"),
        step("1–3  Prepare", "unpack; gaps removed\nfrom every species and\nancestral sequence"),
        step("4  ORF presence", "tBLASTn of each protein vs\nevery species and ancestor;\nintact ORF: Met to stop\n≥ 80% of query, ≥ 40% identity",
             note="→ oldest ORF node"),
        step("5  Synteny", "≥ 40% of the focal\nCDS aligned", note="→ oldest syntenic node"),
    ], [
        step("6–6b  Disablers", "indels and premature\nstops shared by ≥ 2\nnodes between the ORF\nand synteny nodes",
             note="events on ORF-bearing\nnodes dropped"),
        step("7  Classify", "DNG: enabler on the\nbranch before the ORF\nnode, absent from the\nORF clade; else PDNG", note="ORF window set in\ntree_config.json"),
        step("7c–7e  Corroborate", "both enabler nodes seen\nin an extant species,\nand they are the parent\nand grandparent of the\nORF node; else PDNG"),
        step("8–9  Summarize", "counts by node and\nspecies; loci placed on\ninternal nodes by\ncoordinates"),
        io("Output", "DNG and PDNG catalogs;\nnode summaries\nper species\n→ stage 5"),
    ]], None),
    (5, "Births, deaths and rates", [[
        io("Input", "stage 3 coordinates;\nstage 4 summaries and\ntBLASTn; HOG table;\nproteins; branch\nlengths; GFF3 per species"),
        step("1  Coordinates", "each locus in every\ngenome where tBLASTn\nfinds its ORF", note="stops if gene IDs differ"),
        step("2  Ancestral loci", "genes at the same aligned\nposition joined within the\nORF-node clade; oldest\nORF node; DNG if any\nmember is"),
        step("3  Gene overlap", "removes loci whose CDS\noverlaps the CDS of a\nsame-strand gene; antisense,\nintronic and intergenic\nloci kept"),
    ], [
        step("4  Unclustered duplicates", "loci in no orthogroup:\nBLASTp within species,\nE ≤ 10$^{-5}$; the longest\nprotein of a family kept"),
        step("5  Origins", "loci sharing a HOG form one\nDNG and one PDNG origin\nper lineage; PDNG paralogs\nof a DNG dropped"),
        step("5  Births and deaths", "birth on the branch of the\noldest locus; deaths by\nDollo parsimony on\ntBLASTn presence"),
        io("Output", "births, deaths and\nrates (count / branch\nlength) per branch;\nfate of every gene"),
    ]], None),
]


def nlines(s):
    return s.count("\n") + 1 if s else 0


def box_h(n):
    return 2 * PAD + LINE * (1 + nlines(n["body"]) + nlines(n["note"])) + 0.06


def draw(stage, name, rows, footer):
    parallel = stage == 2
    heights = [max(box_h(n) for n in r) for r in rows]
    foot_h = (nlines(footer) * LINE + 0.15) if footer else 0
    H = sum(heights) + GAP_Y * (len(rows) - 1) + foot_h + 0.1
    assert H <= MAX_H, f"stage {stage}: figure height {H:.2f} in exceeds {MAX_H}"
    fig = plt.figure(figsize=(FIG_W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, FIG_W); ax.set_ylim(0, H); ax.axis("off")
    texts, boxes = [], []

    def arrow(p0, p1, **kw):
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=16, lw=1.5, color=HUE, zorder=3,
                                     shrinkA=0, shrinkB=0, **kw))

    def line(xs, ys):
        ax.plot(xs, ys, lw=1.5, color=HUE, zorder=3, solid_capstyle="butt")

    def box(n, x, ytop, w, h):
        fill = IOFILL if n["kind"] == "io" else "white"
        ls = (0, (4, 2.5)) if n["kind"] == "opt" else "-"
        p = FancyBboxPatch((x, ytop - h), w, h, boxstyle="round,pad=0,rounding_size=0.10",
                           facecolor=fill, edgecolor=HUE if n["kind"] != "io" else NOTE, lw=1.6, linestyle=ls, zorder=2)
        ax.add_patch(p)
        boxes.append((x, ytop - h, w, h))
        y = ytop - PAD - LINE / 2
        t = ax.text(x + w / 2, y, n["title"], ha="center", va="center", fontsize=FS_T, fontweight="bold",
                    color=HUE if n["kind"] != "io" else INK, zorder=5)
        texts.append((t, len(boxes) - 1))
        y -= LINE * 0.5
        if n["body"]:
            t = ax.text(x + w / 2, y - 0.03, n["body"], ha="center", va="top", fontsize=FS_B, color=INK, zorder=5,
                        linespacing=1.22)
            texts.append((t, len(boxes) - 1))
            y -= LINE * nlines(n["body"])
        if n["note"]:
            t = ax.text(x + w / 2, y - 0.05, n["note"], ha="center", va="top", fontsize=FS_B, style="italic",
                        color=NOTE, zorder=5)
            texts.append((t, len(boxes) - 1))

    ytop = H - 0.05
    geo = []
    for r, h in zip(rows, heights):
        k = len(r)
        if parallel and k == 1:
            w = 7.0
            xs = [(FIG_W - w) / 2]
        else:
            w = (FIG_W - 0.1 - GAP_X * (k - 1)) / k
            xs = [0.05 + i * (w + GAP_X) for i in range(k)]
        for n, x in zip(r, xs):
            box(n, x, ytop, w, h)
        geo.append((xs, w, ytop, h))
        if not (parallel and k != 4):
            for i in range(k - 1):
                y = ytop - h / 2
                arrow((xs[i] + w, y), (xs[i + 1], y))
        ytop -= h + GAP_Y

    for a, b in zip(geo, geo[1:]):
        (xa, wa, ta, ha), (xb, wb, tb, hb) = a, b
        ybot, ynext = ta - ha, tb
        ymid = (ybot + ynext) / 2
        if parallel and len(xa) == 1:
            cx = xa[0] + wa / 2
            line([cx, cx], [ybot, ymid])
            line([xb[0] + wb / 2, xb[-1] + wb / 2], [ymid, ymid])
            for x in xb:
                arrow((x + wb / 2, ymid), (x + wb / 2, ynext))
        elif parallel:
            for x in xa:
                line([x + wa / 2, x + wa / 2], [ybot, ymid])
            line([xa[0] + wa / 2, xa[-1] + wa / 2], [ymid, ymid])
            sx = xb[0] + wb / 2
            arrow((sx, ymid), (sx, ynext))
        else:
            ex = xa[-1] + wa / 2
            sx = xb[0] + wb / 2
            line([ex, ex], [ybot, ymid])
            line([sx, ex], [ymid, ymid])
            arrow((sx, ymid), (sx, ynext))

    if footer:
        t = ax.text(FIG_W / 2, 0.08, footer, ha="center", va="bottom", fontsize=FS_B, color=NOTE, linespacing=1.22)
        texts.append((t, None))

    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    bad = []
    for t, bi in texts:
        assert t.get_fontsize() >= MIN_FS
        bb = inv.transform(t.get_window_extent(rend))
        (x0, y0), (x1, y1) = bb
        if bi is None:
            assert x0 >= 0 and x1 <= FIG_W and y0 >= 0, f"stage {stage}: footer outside the figure"
            continue
        bx, by, bw, bh = boxes[bi]
        m = 0.04
        if not (x0 >= bx + m and x1 <= bx + bw - m and y0 >= by + m * 0.5 and y1 <= by + bh - m * 0.5):
            bad.append(f"stage {stage}: {t.get_text()[:40]!r} needs {x1 - x0:.2f} in, box {bw - 2 * m:.2f} in")
    assert not bad, "\n".join(bad)
    out = os.path.join(OUTDIR, f"stage{stage}_workflow")
    fig.savefig(out + ".png", dpi=300)
    fig.savefig(out + ".svg")
    plt.close(fig)
    return out + ".png", H


def deck(made):
    from pptx import Presentation
    from pptx.oxml.ns import qn
    from pptx.util import Inches, Pt
    from lxml import etree

    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
    blank = prs.slide_layouts[6]
    for stage, name, png, H in made:
        s = prs.slides.add_slide(blank)
        tb = s.shapes.add_textbox(Inches(0.2), Inches(0.12), Inches(12.9), Inches(0.6))
        p = tb.text_frame.paragraphs[0]
        r = p.add_run()
        r.text = f"LINGUA stage {stage} · {name}"
        r.font.size, r.font.bold, r.font.name = Pt(28), True, "Arial"
        top = 0.85 + (MAX_H - H) / 2
        s.shapes.add_picture(png, Inches((13.333 - FIG_W) / 2), Inches(top), width=Inches(FIG_W))
        nb = s.shapes.add_textbox(Inches(12.3), Inches(7.05), Inches(0.9), Inches(0.4))
        para = nb.text_frame.paragraphs[0]
        para.alignment = 3
        fld = etree.SubElement(para._p, qn("a:fld"), id="{B6F15528-21DE-4FAA-801E-634DDDAF4B2B}", type="slidenum")
        rpr = etree.SubElement(fld, qn("a:rPr"), lang="en-US", sz="1400")
        etree.SubElement(fld, qn("a:t")).text = str(stage + 1)
    prs.core_properties.author = prs.core_properties.last_modified_by = "Casola, Claudio"
    prs.core_properties.title = "LINGUA workflows"
    path = os.path.join(OUTDIR, "LINGUA_workflows.pptx")
    prs.save(path)
    return path


def main():
    if not os.path.isdir(os.path.join(os.getcwd(), "5_birth_death")):
        sys.exit("ABORT, run from the LINGUA root")
    os.makedirs(OUTDIR, exist_ok=True)
    targets = [os.path.join(OUTDIR, f"stage{s}_workflow.{e}") for s, *_ in STAGES for e in ("png", "svg")]
    targets.append(os.path.join(OUTDIR, "LINGUA_workflows.pptx"))
    existing = [t for t in targets if os.path.exists(t)]
    if existing and not OVERWRITE:
        sys.exit(f"ABORT, {existing[0]} exists (use --overwrite)")
    made = []
    for stage, name, rows, footer in STAGES:
        png, H = draw(stage, name, rows, footer)
        made.append((stage, name, png, H))
        print(f"stage {stage}: {os.path.relpath(png)} ({H:.2f} in tall)", file=sys.stderr)
    print(f"wrote {os.path.relpath(deck(made))}", file=sys.stderr)


if __name__ == "__main__":
    main()
