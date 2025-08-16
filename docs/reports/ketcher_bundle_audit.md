# Ketcher Bundle Audit

**Search Result:**  
No local Ketcher bundle (`index.html`, CSS, JS, icons) was found in the current codebase using direct path and filename search.

---

## Audit Summary

- **Found?** No complete self-contained Ketcher build detected.
- **Paths checked:** Any directory or file containing "ketcher" and matching index.html, CSS, JS, or icon naming.

---

## Proposal: Minimal Vendor Folder Structure (`iam2.0/static/ketcher`)

If integration of a local (offline) Ketcher build is required, place the full bundle under:

```
iam2.0/
└── static/
    └── ketcher/
        ├── index.html
        ├── css/
        │   └── ketcher.min.css
        ├── js/
        │   ├── ketcher.min.js
        │   └── ketcher.min.js.map
        ├── img/
        │   ├── <all icons: .svg, .png>
        ├── fonts/
        │   └── <font files>
        └── LICENSE
```

### Required Files for a Minimal Self-Contained Ketcher Vendor Bundle

- `index.html` (Ketcher app entrypoint)
- `css/ketcher.min.css` (core styles)
- `js/ketcher.min.js` and optionally `ketcher.min.js.map` (core app JS)
- `img/` (all icons and graphics used by the editor, typically dozens of SVG/PNG files)
- `fonts/` (any custom fonts used by the editor)
- `LICENSE` (open source license for Ketcher)
- Additional assets referenced in HTML/JS (e.g., manifest, favicon)

---

**Action Required:**  
If you wish to use a local Ketcher build, download the official release from [https://github.com/epam/ketcher](https://github.com/epam/ketcher) and place all required files into `iam2.0/static/ketcher` as above.

**Note:**  
If you have a partial bundle, ensure all assets referenced by `index.html` are present—missing JS, CSS, or icons will break the editor.

--- 

**No local bundle was found; full vendorization is recommended for offline use.**
