"""Publication vector schematic for the paper's problem statement.

Plot-only authoring: no FE/mesh calculations. Dimensions and curve shape are
illustrative; no physical values are printed. Global axes have the origin at
mid-height on the left boundary to match the manuscript's physical y frame.

Usage: python generate_specimen_geometry.py [output.pdf]
Requires: reportlab
"""
from pathlib import Path
import sys
from math import cos, sin, pi
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor


OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name("specimen_geometry.pdf")
OUT.parent.mkdir(parents=True, exist_ok=True)
PAGE_W, PAGE_H = 620, 430
C = canvas.Canvas(str(OUT), pagesize=(PAGE_W, PAGE_H), pageCompression=1)
C.setTitle("Plate with an eccentric circular hole: geometry and loading")
C.setAuthor("CrackFromHoleLEFM manuscript")
BLACK = HexColor("#111111")
LIGHT = HexColor("#777777")
C.setStrokeColor(BLACK)
C.setFillColor(BLACK)

# Rectangular plate in the schematic, with the correct 300:200 shape ratio.
x0, y0 = 91.0, 90.0
W, H = 380.0, 253.333333333333
xm, ym = x0 + W / 2, y0 + H / 2
x1, y1 = x0 + W, y0 + H
sc = W / 300.0
xc = x0 + 170 * sc
yc = ym - 20 * sc
r = 30 * sc

C.setLineWidth(1.3)
C.rect(x0, y0, W, H)
C.circle(xc, yc, r, stroke=1, fill=0)

# Thin center locator from the global origin and vertical reference level.
C.setDash(3, 3)
C.setLineWidth(0.65)
C.setStrokeColor(LIGHT)
C.line(x0, yc, xc, yc)
C.line(xc, yc, xc, ym)
C.setDash()
C.setStrokeColor(BLACK)

# Global coordinate axes (left side midpoint, original y=0).
originx, originy = x0, ym
C.setLineWidth(0.8)
C.line(originx, originy, originx + 42, originy)
C.line(originx, originy, originx, originy + 43)

def arrow_line(xa, ya, xb, yb, head=6, span=3.1, lw=0.9):
    import math
    ang = math.atan2(yb - ya, xb - xa)
    C.setLineWidth(lw)
    C.line(xa, ya, xb, yb)
    for s in (-1,1):
        px = xb - head*math.cos(ang) + s*span*math.sin(ang)
        py = yb - head*math.sin(ang) - s*span*math.cos(ang)
        C.line(xb, yb, px, py)

def double_arrow(xa,ya,xb,yb,head=5.2,span=2.3,lw=0.7):
    arrow_line(xa,ya,xb,yb,head,span,lw)
    arrow_line(xb,yb,xa,ya,head,span,lw)

# Axis arrowheads: previous thin line first, then tips.
arrow_line(originx + 26, originy, originx + 45, originy, head=6, lw=0.8)
arrow_line(originx, originy + 25, originx, originy + 47, head=6, lw=0.8)
C.setFont("Times-Italic", 13)
C.drawString(originx + 48, originy - 4, "x")
C.drawString(originx - 5, originy + 54, "y")
C.setFont("Times-Roman", 9.5)
C.drawRightString(originx - 6, originy - 12, "O")

# Hole center and its (x_c,y_c) locator.
C.setFillColor(BLACK)
C.circle(xc, yc, 1.7, stroke=0, fill=1)
C.setLineWidth(0.7)
C.line(xc-1,yc+1,xc-35,yc+37)
C.setFont("Times-Italic", 12)
text_x, text_y = xc-90, yc+43
C.drawString(text_x, text_y, "(x")
C.setFont("Times-Italic", 8.5); C.drawString(text_x+11.0, text_y-3, "c")
C.setFont("Times-Italic", 12); C.drawString(text_x+15.7,text_y,", y")
C.setFont("Times-Italic", 8.5); C.drawString(text_x+28.7,text_y-3,"c")
C.setFont("Times-Italic", 12); C.drawString(text_x+33.7,text_y,")")

# Radius (toward upper right). Endpoint is deliberately away from crack mouth.
ang = pi/3.35
xe,ye=xc+r*cos(ang),yc+r*sin(ang)
arrow_line(xc,yc,xe,ye,head=5,span=2.1,lw=0.78)
C.setFont("Times-Italic", 12)
C.drawString(xc + 22, yc + 13, "R")

# Gently deflected and then flattening crack, purely schematic, not to scale.
phi = -1.5606127781*pi/180
sx, sy = xc+r*cos(phi), yc+r*sin(phi)
C.setLineWidth(1.65)
P=C.beginPath();P.moveTo(sx,sy)
P.curveTo(sx+15,sy-0.3,sx+27,sy-14,sx+58,sy-23)
P.curveTo(sx+78,sy-29,x1-31,sy-26,x1-15,sy-26)
C.drawPath(P)
C.setFont("Times-Italic",11.5)
C.circle(sx,sy,1.6,stroke=0,fill=1)
C.drawString(sx+2,sy+9,"P")
C.setFont("Times-Roman",7.5);C.drawString(sx+9.5,sy+6,"0")
C.setFont("Times-Italic",11.5);C.drawString(x1-24,sy-43,"P")
C.setFont("Times-Italic",8);C.drawString(x1-16,sy-46,"k")

# Uniform outward tensile loading; same unsigned sigma on both boundaries.
for i in range(13):
    x=x0+12 + (W-24)*i/12
    arrow_line(x,y1,x,y1+29,head=5,span=2.3,lw=0.72)
    arrow_line(x,y0,x,y0-29,head=5,span=2.3,lw=0.72)
# Unlabelled tensile magnitude sigma: small vector glyph to avoid
# dependent/embedded symbol fonts or plus/minus signs.
def draw_sigma(x,y):
    P=C.beginPath()
    P.moveTo(x+13,y+11)
    P.curveTo(x+9,y+13,x+4,y+12,x+2,y+9)
    P.curveTo(x-1,y+5,x+1,y+1,x+5,y+1)
    P.curveTo(x+10,y+1,x+12,y+7,x+10,y+10)
    C.setLineWidth(1.05)
    C.drawPath(P)
draw_sigma(x1+11,y1+17)
draw_sigma(x1+11,y0-30)

# Overall width 2A (below the lower row of traction arrows).
C.setStrokeColor(LIGHT);C.setLineWidth(.65)
C.line(x0,y0-7,x0,17)
C.line(x1,y0-7,x1,17)
C.setStrokeColor(BLACK)
double_arrow(x0,21,x1,21)
C.setFillColor(HexColor("#ffffff"));C.rect(xm-18,13,36,17,fill=1,stroke=0)
C.setFillColor(BLACK);C.setFont("Times-Italic",14);C.drawCentredString(xm,17,"2A")

# Overall height B (on the right).
C.setStrokeColor(LIGHT);C.line(x1+7,y0,x1+42,y0);C.line(x1+7,y1,x1+42,y1)
C.setStrokeColor(BLACK)
double_arrow(x1+37,y0,x1+37,y1)
C.setFillColor(HexColor("#ffffff"));C.rect(x1+28,ym-10,19,20,fill=1,stroke=0)
C.setFillColor(BLACK);C.setFont("Times-Italic",14);C.drawCentredString(x1+37,ym-5,"B")

C.showPage();C.save()
print(str(OUT))
