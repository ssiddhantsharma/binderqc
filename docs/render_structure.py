"""Render a binder-target complex as a cartoon under a translucent molecular
surface (binder blue, target grey), via 3Dmol.js in headless Chromium. Real
coordinates, not a drawing. Works on any PDB/CIF, so it is handy for figures
beyond this repo.

    python docs/render_structure.py [structure] [binder_chain] [target_chain] [out.png]

Defaults render the bundled 7JZU complex (chain A binder, B target) to docs/structure.png.
Needs: pip install binderqc[docs]  &&  python -m playwright install chromium
"""
import base64
import json
import os
import sys

from playwright.sync_api import sync_playwright

HERE = os.path.dirname(os.path.abspath(__file__))


def render(structure_path, binder_chain="A", target_chain="B", out_path=None, size=1100):
    out_path = out_path or os.path.join(HERE, "structure.png")
    with open(structure_path) as fh:
        text = fh.read()
    fmt = "cif" if structure_path.lower().endswith((".cif", ".mmcif")) else "pdb"
    html = """<!doctype html><html><head>
<script src="https://cdn.jsdelivr.net/npm/3dmol@2.4.2/build/3Dmol-min.js"></script>
</head><body style="margin:0;background:white">
<div id="v" style="width:%dpx;height:%dpx;position:relative"></div>
<script>
const data = %s, fmt = %s, B = %s, T = %s;
(async () => {
 try {
  const v = $3Dmol.createViewer("v", {backgroundColor:"white"});
  v.addModel(data, fmt);
  v.setStyle({}, {});
  v.setStyle({chain:T}, {cartoon:{color:"0x9aa0a6"}});
  v.setStyle({chain:B}, {cartoon:{color:"0x2f6fb0"}});
  v.setViewStyle({style:"outline", color:"0x222222", width:0.04});
  await v.addSurface($3Dmol.SurfaceType.SES, {opacity:0.5, color:"0xbcc2c6"}, {chain:T});
  await v.addSurface($3Dmol.SurfaceType.SES, {opacity:0.4, color:"0x8fb9ea"}, {chain:B});
  v.zoomTo(); v.rotate(90, "x"); v.rotate(15, "y"); v.render();
  await new Promise(r => setTimeout(r, 900));
  window.__png = v.pngURI(); window.__ready = true;
 } catch (e) { window.__err = e.message; window.__ready = true; }
})();
</script></body></html>""" % (size, size, json.dumps(text), json.dumps(fmt),
                              json.dumps(binder_chain), json.dumps(target_chain))
    with sync_playwright() as p:
        b = p.chromium.launch(args=["--use-gl=angle", "--use-angle=swiftshader", "--ignore-gpu-blocklist"])
        pg = b.new_page(viewport={"width": size, "height": size}, device_scale_factor=2)
        pg.set_content(html, wait_until="load")
        pg.wait_for_function("window.__ready === true", timeout=120000)
        err, uri = pg.evaluate("window.__err || ''"), pg.evaluate("window.__png || ''")
        b.close()
    if err:
        raise SystemExit("render failed: " + err)
    open(out_path, "wb").write(base64.b64decode(uri.split(",", 1)[1]))
    print("wrote", out_path)


if __name__ == "__main__":
    a = sys.argv[1:]
    structure = a[0] if len(a) > 0 else os.path.join(HERE, "..", "tests", "data", "7JZU_LCB1_RBD.pdb")
    binder = a[1] if len(a) > 1 else "A"
    target = a[2] if len(a) > 2 else "B"
    out = a[3] if len(a) > 3 else None
    render(structure, binder, target, out)
