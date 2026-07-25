/* ============================================================================
   COSE book-nav (v0.1) — integrates the cross-site project map into the NATIVE
   sidebar of framework docs sites (Jupyter Book / sphinx-book-theme / pydata),
   which already have their own left nav. Additive only: it appends a "CoSE
   Project Family" section from window.BARKER_SITES (sites.js) to the primary
   sidebar. No DOM restructuring, so it can't fight the theme's layout. If the
   sidebar isn't found it does nothing.

   Use: load sites.js then this file via the book's html_js_files.
   ========================================================================== */
(function(){
  "use strict";
  function build(){
    var reg = window.BARKER_SITES;
    if(!reg || !reg.groups || document.getElementById("cose-booknav")) return;
    var target = document.querySelector(".bd-sidebar-primary .sidebar-primary-items__end")   // jupyter book / pydata
      || document.querySelector(".bd-sidebar-primary nav.bd-links")
      || document.querySelector(".md-sidebar--primary .md-nav--primary")                     // mkdocs material
      || document.querySelector(".md-sidebar--primary .md-sidebar__inner")
      || document.querySelector(".bd-sidebar-primary")
      || document.querySelector("nav.bd-links");
    if(!target) return;

    var st = document.createElement("style");
    st.textContent =
      "#cose-booknav{padding:14px 0 8px;border-top:1px solid var(--pst-color-border,#e5e9f0);margin-top:14px}"+
      "#cose-booknav .cn-title{font-size:.72rem;font-weight:700;letter-spacing:.06em;text-transform:uppercase;"+
        "color:var(--pst-color-primary,#3B6EA5);margin:0 0 6px;padding:0 6px}"+
      "#cose-booknav .cn-grp{font-size:.64rem;font-weight:700;letter-spacing:.05em;text-transform:uppercase;"+
        "color:var(--pst-color-text-muted,#8a8a8a);margin:11px 6px 2px}"+
      "#cose-booknav a{display:block;padding:4px 8px;font-size:.82rem;line-height:1.3;border-radius:6px;"+
        "color:var(--pst-color-text-base,inherit);text-decoration:none}"+
      "#cose-booknav a:hover{background:var(--pst-color-surface,rgba(59,110,165,.08));text-decoration:none}";
    document.head.appendChild(st);

    var box = document.createElement("div");
    box.id = "cose-booknav"; box.className = "sidebar-primary-item";
    var t = document.createElement("p"); t.className = "cn-title";
    t.textContent = "CoSE Project Family"; box.appendChild(t);
    reg.groups.forEach(function(g){
      if(!g.items.some(function(it){ return it.url && it.live !== false; })) return;
      var gh = document.createElement("div"); gh.className = "cn-grp";
      gh.textContent = g.name; box.appendChild(gh);
      g.items.forEach(function(it){
        if(!it.url || it.live === false) return;
        var a = document.createElement("a"); a.href = it.url;
        a.textContent = (it.emoji ? it.emoji + " " : "") + it.title;
        box.appendChild(a);
      });
    });
    target.appendChild(box);
  }
  if(document.readyState === "loading") document.addEventListener("DOMContentLoaded", build);
  else build();
  // MkDocs Material (and other SPA-style docs) re-render the sidebar on
  // in-page navigation; re-inject if our section gets removed. build() is a
  // cheap no-op when the section is already present.
  if("MutationObserver" in window){
    new MutationObserver(function(){
      if(!document.getElementById("cose-booknav")) build();
    }).observe(document.body, {childList:true, subtree:true});
  }
})();
