/* ============================================================================
   Barker Lab shared theme — behaviour (v0.3)
   - Top-left COSE brand logo + a "back to hub" link. Cross-site navigation
     lives on the hub pages (personal site → hub → project); the old
     slide-out map rail was retired 2026-10-21.
   - The hub link returns to the hub the visitor came through: a sub-hub links
     its cards with ?cose_scope=NAME, which is remembered per-tab in
     sessionStorage; otherwise it goes to the main CoSE hub.
   - Light/dark toggle (persisted) + COSE link in the footer.
   - Progressive enhancement: if this script doesn't run, the page is unchanged.

   To use on a page:  set  <body data-site-id="THIS_REPO_SLUG">
   and include sites.js, then theme.js, at the end of <body>.
   ========================================================================== */
(function(){
  "use strict";
  var LS_THEME = "barker.theme";
  var slug = document.body.getAttribute("data-site-id") || "";
  var LOGO = document.body.getAttribute("data-brand-logo") || "assets/cose-logo.png";
  var BRAND_URL = document.body.getAttribute("data-brand-url") || "https://cosecloud.com/";
  var reg = window.BARKER_SITES || {};

  // Forget the retired rail's saved state.
  try{ ["barker.map.open","barker.map.panel","barker.rail"].forEach(function(k){ localStorage.removeItem(k); }); }catch(e){}

  /* ---------- which hub did the visitor come from? ---------- */
  var scopes = reg.scopes || {};
  var SCOPE = "";
  try {
    if(slug === "cose-hub"){ sessionStorage.removeItem("cose.scope"); }   // the main hub resets it
    else if(scopes[slug]){ sessionStorage.setItem("cose.scope", slug); }  // a sub-hub sets its own
    var sc = new URLSearchParams(location.search).get("cose_scope");
    if(sc !== null){ if(sc){ sessionStorage.setItem("cose.scope", sc); } else { sessionStorage.removeItem("cose.scope"); } }
    else { sc = sessionStorage.getItem("cose.scope"); }
    if(sc && scopes[sc] && scopes[sc].hub){ SCOPE = sc; }
  } catch(e){}

  var hubUrl = "", hubLabel = "", hubTitle = "";
  if(slug === "cose-hub"){ /* already on the main hub */ }
  else if(scopes[slug]){ hubUrl = reg.hub; hubLabel = "All projects"; hubTitle = "All projects — CoSE hub"; }
  else if(SCOPE){ hubUrl = scopes[SCOPE].hub; hubLabel = scopes[SCOPE].hubLabel || scopes[SCOPE].label;
    hubTitle = "Back to the " + hubLabel + " hub"; }
  else { hubUrl = reg.hub; hubLabel = "CoSE hub"; hubTitle = "Back to the CoSE hub"; }

  /* ---------- top-left bar: COSE brand logo + back-to-hub link ---------- */
  var topbar = el("div", {class:"topbar"});
  var brand = el("a", {class:"cose-brand", href:BRAND_URL, target:"_blank",
    rel:"noopener", title:"COSE — cosecloud.com", "aria-label":"COSE — cosecloud.com"});
  brand.appendChild(el("img", {src:LOGO, alt:"COSE"}));
  topbar.appendChild(brand);
  if(hubUrl){
    var hubLink = el("a", {class:"cose-hub-btn", href:hubUrl, title:hubTitle, "aria-label":hubTitle});
    hubLink.innerHTML = '<svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M19 12H5M11 18l-6-6 6-6"/></svg>'
      + '<span class="lbl-long">' + esc(hubLabel) + '</span><span class="lbl-short">Hub</span>';
    topbar.appendChild(hubLink);
  }

  var SUN = '<svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><circle cx="12" cy="12" r="4"/><path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4"/></svg>';
  var MOON = '<svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 12.8A9 9 0 1 1 11.2 3a7 7 0 0 0 9.8 9.8z"/></svg>';
  var themeBtn = el("button", {class:"cose-theme-btn",
    "aria-label":"Toggle light or dark theme", title:"Toggle light / dark"});
  // In the top bar beside the hub link (theme.css pins it top-right instead).
  if(document.body.getAttribute("data-cose-themetoggle") !== "off"){ topbar.appendChild(themeBtn); }

  /* Wrap existing body content in div.page. The wrapper was added for the
     retired rail; it's kept so the DOM (and any CSS written against it) on
     the ~70 pages using this kit stays exactly as it was. */
  var page = el("div", {class:"page"});
  while(document.body.firstChild){ page.appendChild(document.body.firstChild); }
  document.body.appendChild(topbar);
  document.body.appendChild(page);

  /* COSE brand link in the footer, immediately before the author name.
     Falls back through common footer shapes (radiation uses footer>.wrap>p;
     the shared template uses a bare <footer class="site">text). */
  var footP = page.querySelector("footer .wrap p") || page.querySelector("footer .wrap")
    || page.querySelector("footer p") || page.querySelector("footer");
  if(footP){
    var fLink = el("a", {class:"cose-foot", href:BRAND_URL, target:"_blank",
      rel:"noopener", title:"COSE — cosecloud.com", "aria-label":"COSE — cosecloud.com"});
    fLink.appendChild(el("img", {src:LOGO, alt:"COSE"}));
    fLink.appendChild(document.createTextNode("COSE"));
    footP.insertBefore(fLink, footP.firstChild);
  }

  /* Match the top bar to the host page's ACTUAL background brightness,
     so a site that hard-codes light or dark (ignoring the OS preference) still
     gets a matching bar. Scoped tokens only recolour the CoSE chrome. */
  (function(){
    function lum(c){
      var m = c && c.match(/rgba?\(([^)]+)\)/); if(!m) return null;
      var p = m[1].split(",").map(parseFloat);
      if(p.length >= 4 && p[3] === 0) return null;            // transparent
      return (0.2126*p[0] + 0.7152*p[1] + 0.0722*p[2]) / 255;
    }
    var l = lum(getComputedStyle(document.body).backgroundColor);
    if(l === null) l = lum(getComputedStyle(document.documentElement).backgroundColor);
    if(l === null) l = window.matchMedia("(prefers-color-scheme: dark)").matches ? 0 : 1;
    topbar.classList.add(l < 0.45 ? "cose-dark" : "cose-light");
  })();

  /* ---------- light / dark theme toggle ----------
     Sets data-theme on <html>; cose-theme.css / theme.css define the token
     values for html[data-theme="light|dark"], overriding the page's own
     prefers-color-scheme media queries. Persisted; default follows the OS. */
  function effectiveTheme(){
    return document.documentElement.getAttribute("data-theme")
      || (window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
  }
  function applyTheme(t){
    document.documentElement.setAttribute("data-theme", t);
    themeBtn.innerHTML = t === "dark" ? SUN : MOON;
    themeBtn.setAttribute("title", t === "dark" ? "Switch to light" : "Switch to dark");
    topbar.classList.remove("cose-dark","cose-light");
    topbar.classList.add(t === "dark" ? "cose-dark" : "cose-light");
    try{ localStorage.setItem(LS_THEME, t); }catch(e){}
  }
  themeBtn.addEventListener("click", function(){
    applyTheme(effectiveTheme() === "dark" ? "light" : "dark"); });
  var savedTheme = null;
  try{ savedTheme = localStorage.getItem(LS_THEME); }catch(e){}
  if(savedTheme){ applyTheme(savedTheme); }
  else { themeBtn.innerHTML = effectiveTheme() === "dark" ? SUN : MOON; }

  /* ---------- tiny helpers ---------- */
  function el(tag, attrs){ var n=document.createElement(tag);
    if(attrs) for(var k in attrs){ n.setAttribute(k, attrs[k]); }
    return n; }
  function esc(s){ return String(s).replace(/[&<>"]/g,function(c){
    return {"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]; }); }
})();
