/* ============================================================================
   COSE — shared site registry (the "site map" content + future hub source).
   Edit THIS one file to add/rename a project; every page that ships a copy
   (or loads the hosted copy) gets the updated cross-site nav.

   Groups + titles follow Richard's authoritative list (2026-07-24).
   `id` must match the repo slug (last path segment of the github.io URL) so a
   page can mark "you are here" via <body data-site-id="…">.
   `live:false` renders as "page pending" and is skipped by the hub.
   ========================================================================== */
window.BARKER_SITES = {
  brand: { name:"COSE", url:"https://cosecloud.com/", logo:"assets/cose-logo.png" },
  hub: "https://dr-richard-barker.github.io/CoSE_Cloud/Hub/",
  groups: [
    {
      name: "Featured",
      items: [
        { id:"LunarLeaf-CFD", emoji:"🍃", title:"Lunar LEAF — Photorespiration",
          desc:"CFD model of photorespiration in a lunar growth chamber",
          url:"https://dr-richard-barker.github.io/LunarLeaf-CFD/" },
        { id:"deepspace-seed-stress-decoder", emoji:"🌰", title:"DeepSpace Seed Stress Decoder",
          desc:"Decoding seed stress signatures for deep-space conditions",
          url:"https://dr-richard-barker.github.io/deepspace-seed-stress-decoder/" },
        { id:"Plant_response_to_radiation", emoji:"☢️", title:"OSDR Radiation Review",
          desc:"Radiation review & kinetic pattern-recognition across NASA OSDR",
          url:"https://dr-richard-barker.github.io/Plant_response_to_radiation/" },
        { id:"Astronaut_health_search", emoji:"🧬", title:"Astronaut Health Summary",
          desc:"Searchable summary of astronaut health data",
          url:"https://dr-richard-barker.github.io/Astronaut_health_search/" },
        { id:"Astronaut_trends", emoji:"📊", title:"Astronaut Trends Dashboard",
          desc:"Interactive dashboard of astronaut health trends",
          url:"https://dr-richard-barker.github.io/Astronaut_trends/" },
      ]
    },
    {
      name: "Spaceflight Omics & Transcriptomics",
      items: [
        { id:"Tropism_autodecoder_2026", emoji:"🧭", title:"Tropism Autodecoder 2026",
          desc:"Auto-decoder atlas of plant tropism responses",
          url:"https://dr-richard-barker.github.io/Tropism_autodecoder_2026/" },
        { id:"arabidopsis-spaceflight-omics", emoji:"🌿", title:"Arabidopsis Spaceflight Omics",
          desc:"LASSO biomarkers, single-cell atlas integration & Ca²⁺ cell–cell signalling across NASA OSDR",
          url:"https://dr-richard-barker.github.io/arabidopsis-spaceflight-omics/" },
        { id:"B_rappa_LLGCSS", emoji:"🌸", title:"B. rapa — Floral Scent × Radiation",
          desc:"Does galactic cosmic radiation alter floral scent? WIP transcriptomics in Brassica rapa",
          url:"https://dr-richard-barker.github.io/B_rappa_LLGCSS/" },
        { id:"APEX05_results_and_code", emoji:"🪴", title:"APEX-05 Clean-Up & Analysis",
          desc:"APEX-05 results and reproducible analysis code",
          url:"https://dr-richard-barker.github.io/APEX05_results_and_code/" },
        { id:"TICTOC", emoji:"🧵", title:"TICTOC Project Clean-Up",
          desc:"TICTOC project data clean-up and documentation",
          url:"https://dr-richard-barker.github.io/TICTOC/" },
        { id:"smallRNAseq-DREAM", emoji:"🧫", title:"MicroRNA Analysis Pipeline",
          desc:"Cross-species small-RNA-seq (miRNA) pipeline & OSDR mining test",
          url:"https://dr-richard-barker.github.io/smallRNAseq-DREAM/" },
      ]
    },
    {
      name: "Microbiome & Multi-Omics Reviews",
      items: [
        { id:"osdr-plant-microbiome", emoji:"🦠", title:"OSDR Plant Microbiome Review",
          desc:"Plant microbiome review & manuscript",
          url:"https://dr-richard-barker.github.io/osdr-plant-microbiome/" },
        { id:"veg05-integrated-omics", emoji:"🍅", title:"VEG-05 Integrated-Omics",
          desc:"Multi-omics of ISS dwarf tomato (VEG-05): red- vs blue-rich lighting vs KSC controls",
          url:"https://dr-richard-barker.github.io/veg05-integrated-omics/" },
        { id:"aph-physiospace", emoji:"🧠", title:"APH Tissue-Specific PhysioSpace DL",
          desc:"Tissue-specific PhysioSpace deep learning & salicylic-acid comparison",
          url:"", live:false },
      ]
    },
    {
      name: "Interactive Notebooks, Web Tools & Pipelines",
      items: [
        { id:"Anthocyanin-Image-analysis", emoji:"🍁", title:"Anthocyanin Image Analysis Tool",
          desc:"Browser tool for anthocyanin image analysis",
          url:"https://dr-richard-barker.github.io/Anthocyanin-Image-analysis/" },
        { id:"astroroot", emoji:"🛰️", title:"AstroRoot",
          desc:"In-browser root image analysis for the classroom",
          url:"https://dr-richard-barker.github.io/astroroot/" },
        { id:"virtual-root", emoji:"🌱", title:"Virtual Root",
          desc:"Interactive auxin-transport root model",
          url:"https://dr-richard-barker.github.io/virtual-root/" },
        { id:"madwest-astrobotany", emoji:"🚀", title:"MadWest Astrobotany",
          desc:"MadWest astrobotany outreach & resources",
          url:"https://dr-richard-barker.github.io/madwest-astrobotany/" },
      ]
    },
    {
      name: "Games & 3D Simulations",
      items: [
        { id:"Settlers_of_Mars_3D_LLM", emoji:"🪐", title:"Settlers of Mars (3D + LLM)",
          desc:"Space-biology education game: browser colony sim, AI role-play & a 3D-printed board",
          url:"https://dr-richard-barker.github.io/Settlers_of_Mars_3D_LLM/" },
      ]
    },
    {
      name: "Documentation & Education Hubs",
      items: [
        { id:"OSDR_jupyter_book.io", emoji:"📓", title:"OSDR Jupyter Book (TOAST10)",
          desc:"OSDR Jupyter Book / TOAST10 interactive notebooks",
          url:"https://dr-richard-barker.github.io/OSDR_jupyter_book.io/" },
        { id:"AIRI", emoji:"🌍", title:"Astrobotany International Research Initiative",
          desc:"AIRI — astrobotany international research initiative",
          url:"https://dr-richard-barker.github.io/AIRI/" },
      ]
    }
  ]
};
