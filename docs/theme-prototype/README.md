# Barker Lab — uniform GitHub Pages theme (plan + working prototype)

**Goal:** give every one of Richard's project GitHub Pages the same minimalist look as
[Plant_response_to_radiation](https://dr-richard-barker.github.io/Plant_response_to_radiation/),
and add an **optional, toggleable left-hand "map"** that combines a **document map** (this page's
sections) and a **site map** (links to all the other project pages).

This folder is a **non-destructive prototype**: the live page at `docs/index.html` is untouched.
The prototype lives at `docs/theme-prototype/` so you can preview it side-by-side before we commit
to a rollout.

> **Status:** prototype built and verified in-browser on 2026-07-24. The cross-site registry
> (`sites.js`) is a *seed list* built from repos found on this machine — **titles, URLs, and
> which pages are actually live still need Richard's confirmation** before rollout.

---

## 1. What we're keeping (the reference design)

The radiation page is a single self-contained `index.html` — **no Jekyll, no build step**
(`.nojekyll` is present). Its design system is worth preserving verbatim:

| Element | Value |
|---|---|
| Content column | `max-width: 960px`, centered, 20px gutters |
| Colour tokens | CSS custom properties (`--bg --fg --muted --line --card --accent --accent2`) |
| Light/dark | automatic via `@media (prefers-color-scheme: dark)` |
| Accent | `#3B6EA5` blue + `#3FB6A8` teal |
| Type | system font stack, 16px/1.6 |
| Components | eyebrow, hero, `.btn`, `.cards`/`.card`/`.stat`, `figure`, `table`, `.disc` callout |

The prototype lifts these tokens and component styles **unchanged** into a shared `theme.css`, so a
themed page is pixel-identical to the current one when the map is closed.

## 2. What the prototype adds

A portable, dependency-free **theme kit** — four files, no framework, no build:

```
docs/theme-prototype/
├── index.html              # the radiation page, re-pointed at the shared theme
└── assets/
    ├── theme.css           # design tokens + components (verbatim) + left-rail + COSE brand styles
    ├── theme.js            # builds the rail, auto document map, scroll-spy, toggle, COSE branding, persistence
    ├── sites.js            # the SHARED cross-site registry (the "site map" content)
    └── cose-logo.png       # COSE brand mark (from cosecloud.com), 256×256
```

**Behaviour**
- A small **"Map" toggle** sits top-left. Closed by default → page looks exactly as it does today.
- Open reveals a left rail with two tabs:
  - **On this page** — a document map auto-generated from each `<section>`'s `<h2>`
    (no manual list to maintain), with **scroll-spy** highlighting the current section.
  - **All projects** — the site map, grouped by theme, with the current site highlighted.
- **Responsive:** on wide screens (≥1180px) the open rail *pushes* the content; on narrow screens it
  *overlays* with a dim scrim (tap-to-close, Esc-to-close).
- **Remembers** open/closed state and last-used tab in `localStorage`.
- **COSE branding:** the COSE logo appears **top-left next to the Map button** and **in the footer
  before the author name** — both link to <https://cosecloud.com/>. Injected automatically by
  `theme.js`; logo path and target are overridable per page via `data-brand-logo` / `data-brand-url`
  on `<body>` (defaults: `assets/cose-logo.png`, `https://cosecloud.com/`).
- **Progressive enhancement:** if `theme.js` fails to load, the page still renders normally.

**Per-page cost to adopt:** link three files and set one attribute —
```html
<link rel="stylesheet" href="assets/theme.css">
...
<body data-site-id="THIS_REPO_SLUG">
...
<script src="assets/sites.js"></script>
<script src="assets/theme.js"></script>
```
The document map needs **nothing else** — it reads the page's own `<section>`/`<h2>` structure.

## 3. Preview it

```bash
cd docs/theme-prototype
python3 -m http.server 8799
# open http://localhost:8799/index.html — toggle the "Map" button top-left
```

## 4. The sites this covers (Richard's authoritative list, 2026-07-24)

17 live sites in 4 groups (see `sites.js`). Classified by cloning each repo and inspecting its
Pages source. **Five different hosting architectures** — the rollout must handle each differently:

| Site (repo) | Architecture | Rollout effort |
|---|---|---|
| `Plant_response_to_radiation` | static `docs/index.html` + `.nojekyll` | ✅ reference (done) |
| `Tropism_autodecoder_2026` | static `docs/index.html` + `.nojekyll` | 🟢 drop-in |
| `Astronaut_health_search` | static `docs/index.html` + `.nojekyll` | 🟢 drop-in |
| `APEX05_results_and_code` | static `docs/index.html` + `.nojekyll` | 🟢 drop-in |
| `veg05-integrated-omics` | static `docs/index.html` + `.nojekyll` | 🟢 drop-in |
| `Astronaut_trends` | static `docs/index.html` | 🟢 drop-in |
| `LunarLeaf-CFD` | static **root** `index.html` | 🟢 drop-in (root paths) |
| `Anthocyanin-Image-analysis` | static **root** `index.html` | 🟢 drop-in (root paths) |
| `madwest-astrobotany` | static **root** `index.html` + `.nojekyll` | 🟢 drop-in (root paths) |
| `deepspace-seed-stress-decoder` | **Jekyll** (`docs/index.md`) | 🟡 edit Jekyll layout |
| `AIRI` | **GitBook** (`docs/index.md`) | 🟠 GitBook has its own theme system |
| `TICTOC` | **GitBook** + a custom `docs/index.html` landing | 🟠 theme the landing; GitBook body separate |
| `smallRNAseq-DREAM` | **built → gh-pages** (no page source on `main`) | 🔴 edit source template + rebuild |
| `osdr-plant-microbiome` | **built → gh-pages** (no page source on `main`) | 🔴 edit source template + rebuild |
| `OSDR_jupyter_book.io` | **Jupyter Book → gh-pages** | 🔴 theme via `_config.yml`/templates + rebuild |
| `OSDR_plant_spaceflight_omics` | ⚠️ repo **"not found"** (private or renamed) | ❓ needs access |
| `aph-physiospace` (APH PhysioSpace DL) | not published yet ("coming soon") | ⬜ scaffold later |

Legend: 🟢 identical drop-in · 🟡 layout edit · 🟠 bespoke (GitBook) · 🔴 build-pipeline edit + rebuild.

## 5. Rollout approach (by architecture)

The theme kit is plain CSS/JS, so it *can* drop into all of these — but the mechanics differ:

1. **Static pages (9 sites)** — the bulk, lowest risk. Add `assets/theme.{css,js}` + `sites.js` +
   `cose-logo.png`, add the `<link>`/`<script>` tags + `data-site-id`. Root-index sites use
   root-relative asset paths. **Do these first, in batches, verifying each.**

2. **Jekyll (deepspace-seed)** — add the kit to the default layout (`_layouts/default.html`) so every
   rendered page emits `<section><h2>` for the doc-map. Content stays in markdown.

3. **GitBook (AIRI, TICTOC)** — GitBook ships its own theme; the CoSE map/brand can't be injected the
   same way. Options: (a) theme only a custom static landing page (TICTOC already has one), or
   (b) add COSE logo + cross-links via GitBook's own customisation. Decide per site.

4. **Built sites (smallRNAseq-DREAM, osdr-plant-microbiome, OSDR_jupyter_book.io)** — the live HTML is
   generated to a `gh-pages` branch from notebooks/markdown. Theming means editing the build config
   (e.g. Jupyter Book `_config.yml`/`_templates`) and re-running the build, not swapping HTML.

5. **Blocked / pending** — `OSDR_plant_spaceflight_omics` clone returns "not found" (confirm it's
   public / the right name); `aph-physiospace` has no page yet (scaffold from the radiation template
   using real content when ready).

### Keeping the site map in sync
`sites.js` is the single source of truth for cross-site nav. Two options — recommend deciding early:
- **(A) Copy per repo** — simplest, but editing the project list means updating N copies.
- **(B) Host one canonical `sites.js`** on one Pages site (e.g. `Space_Biology_Education.io`) and have
  every page load *that* URL. One edit updates every site's map. Slight coupling / cross-origin load.
  **Recommended** once the list stabilises; a small sync script can also push copies as a fallback.

## 6. Suggested sequence

1. ✅ Prototype on the radiation page + COSE branding *(done — this folder)*.
2. ✅ Authoritative `sites.js` (17 sites, 4 groups) *(done)*.
3. ✅ Decisions locked: **direct-to-live** for simple static sites, **preview path** for tricky ones;
   **copy `sites.js` per repo** (sync via script).
4. 🟡 Static sites:
   - ✅ **Live** (direct-to-`main`, verified): `Astronaut_health_search`, `APEX05_results_and_code`,
     `veg05-integrated-omics`, `Tropism_autodecoder_2026`.
   - ✅ **Preview built** (`index.cose.html` beside the live page, live page untouched):
     `madwest-astrobotany` → `/index.cose.html`; `Astronaut_trends` → `/index.cose.html`
     (dashboard; 14-entry doc-map; rail auto-matches its light theme).
   - ⚠️ **`LunarLeaf-CFD` and `Anthocyanin-Image-analysis` are React SPAs** (`<div id="root">`,
     built bundles) — a static overlay can't theme them. They need the CoSE map/brand as a React
     component + a rebuild. Blocked on: how they deploy (gh-pages branch? Actions build?).
   - `theme.js` now samples the host page's background brightness and matches the rail to a
     hard-coded light/dark site (fixes light-rail-on-dark-page). Older shipped copies (the 4 live
     sites) can be re-synced with this build; no visual change since they follow the OS preference.
5. ⬜ Jekyll layout for `deepspace-seed`.
6. ⬜ GitBook sites (`AIRI`, `TICTOC`) — decide per-site.
7. ⬜ Built sites (`smallRNAseq-DREAM`, `osdr-plant-microbiome`, `OSDR_jupyter_book.io`) — build-config edits.
8. ⬜ Resolve `OSDR_plant_spaceflight_omics` access; scaffold `aph-physiospace` when ready.
9. ⬜ **Hub site** linking all of them (built from `sites.js` — same theme, grid of project cards).
10. ⬜ Promote `theme-prototype/` over the radiation `docs/index.html`.

### Theme kit files (copy set for existing styled pages)
`assets/cose-theme.css` (token override + rail + brand overlay) · `assets/theme.js`
(doc-map from `<h2>`, scroll-spy, toggle, brand, persistence) · `assets/sites.js` (shared registry)
· `assets/cose-logo.png`. Wire-up per page: add the `cose-theme.css` `<link>` after the page's own
stylesheet, set `<body data-site-id="REPO_SLUG">`, and add the two `<script>` tags before `</body>`.

## 7. The hub (step 10)

Once the family shares the theme, a new landing site (e.g. a `cose-hub` repo or a page on
`cosecloud.com`) renders the same `sites.js` groups as a card grid — one entry per project, COSE-
branded, same left map. Because it reads the shared registry, adding a project to `sites.js` adds it
to both every site's map *and* the hub with no extra edits.

## 8. Open decisions for Richard

- **Apply style:** preview path first (recommended — matches "test") or straight to each live page?
- **Sync model:** per-repo copy of `sites.js` (simple) vs one hosted canonical copy (one edit updates all).
- **Sidebar default:** closed on first visit (current) or open on desktop?
- **Scope of unification:** just the left map + COSE brand, or also a shared header/footer?
