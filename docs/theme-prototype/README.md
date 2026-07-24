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

## 4. The repos this would cover

Found locally on 2026-07-24. **They are not uniform today** — three different hosting styles,
which is the main thing the rollout has to reconcile:

| Repo | Current page style | Notes |
|---|---|---|
| `Plant_response_to_radiation` | self-contained `docs/index.html` | the reference; already ideal |
| `Tropism_autodecoder_2026` | self-contained `docs/index.html` | closest to reference |
| `tropism-autodecoder-webtool` | self-contained `docs/index.html` | closest to reference |
| `astroroot` | root `index.html` | not in `docs/` |
| `Space_Biology_Education.io` | **Jekyll** (`_config.yml`) | theme via Jekyll layout |
| `bloodbowl` (Training_LLM…) | **Jekyll** (`_config.yml`) + `game/` | game page separate |
| `deepspace-seed-stress-decoder` | **Jekyll** (`docs/_config.yml`, `index.md`) | markdown-driven |
| `B_rappa_LLGCSS` | no page yet | docs/ has manuscript md only |
| `VEGGIE_Tom_Red_Blue…` | no page yet | — |
| `PhysioSpace_stress_decoding_VEG05` (DeepLearning_VEG05) | no page yet | rename pending |
| `smallRNAseq-DREAM` | no page yet | — |

> ⚠️ This list is what's on disk, not confirmed live Pages. **Please confirm the full/authoritative
> list** (you mentioned you'd send it) and which are actually published — I'll reconcile `sites.js`
> against it.

## 5. Rollout approach (by hosting style)

The theme kit is plain CSS/JS, so it drops into all three styles — but the mechanics differ:

1. **Self-contained static pages** (radiation, both tropism pages, astroroot)
   → copy `assets/theme.{css,js}` + `sites.js`, swap the page's `<style>` block for the stylesheet
   link, add the two scripts and `data-site-id`. Lowest risk; do these first.

2. **Jekyll sites** (education, bloodbowl, deepspace-seed)
   → add `theme.css` and the scripts to the Jekyll **default layout** (`_layouts/default.html`) and
   ensure each rendered page emits `<section><h2>…` structure (or point the doc-map at Jekyll's
   heading output). Content stays in markdown; only the layout changes.

3. **Repos with no page yet** (B_rapa, VEGGIE, PhysioSpace, smallRNAseq)
   → scaffold a new `docs/index.html` from the radiation template, populated from each repo's real
   README/manuscript content (no invented content).

### Keeping the site map in sync
`sites.js` is the single source of truth for cross-site nav. Two options — recommend deciding early:
- **(A) Copy per repo** — simplest, but editing the project list means updating N copies.
- **(B) Host one canonical `sites.js`** on one Pages site (e.g. `Space_Biology_Education.io`) and have
  every page load *that* URL. One edit updates every site's map. Slight coupling / cross-origin load.
  **Recommended** once the list stabilises; a small sync script can also push copies as a fallback.

## 6. Suggested sequence

1. ✅ Prototype on the radiation page *(done — this folder)*.
2. ⬜ Richard confirms the authoritative page list + live URLs → finalise `sites.js`.
3. ⬜ Decide sync model (A copy vs B hosted).
4. ⬜ Roll to the 3 other static pages (lowest risk), verify each in-browser.
5. ⬜ Adapt the Jekyll layouts for the 3 Jekyll sites.
6. ⬜ Scaffold pages for the 4 repos that have none.
7. ⬜ Promote `theme-prototype/` to replace `docs/index.html` on the radiation repo.

## 7. Open questions for Richard

- The **full list of pages** you want unified (with live URLs) — several repos on disk may not have
  Pages enabled.
- Sidebar default: **closed** on first visit (current prototype) or **open** on desktop?
- Should the site map group headings (Radiation / Tropism / VEGGIE / Education) match how you think
  about the projects, or do you want a flat alphabetical list?
- Do you want a shared **header/footer** (lab name, links) unified too, or *only* the left map?
