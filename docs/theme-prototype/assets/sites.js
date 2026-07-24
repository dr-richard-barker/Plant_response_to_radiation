/* ============================================================================
   Shared site registry — the "site map" content, identical across all repos.
   This is the ONE file you edit when you add/rename a project; every page that
   ships a copy (or links the hosted copy) gets the updated cross-site nav.

   NOTE: this is a prototype seed list built from the local repos found on
   2026-07-24. URLs/titles/live-status must be confirmed with Richard before
   rollout — several repos may not have Pages enabled yet.
   ========================================================================== */
window.BARKER_SITES = {
  // slug of the current site is matched against `id` to mark "you are here"
  groups: [
    {
      name: "Radiation & stress",
      items: [
        { id:"Plant_response_to_radiation", title:"Plant Response to Radiation",
          desc:"Kinetic transcriptomics of ionizing-radiation response",
          url:"https://dr-richard-barker.github.io/Plant_response_to_radiation/", live:true },
        { id:"deepspace-seed-stress-decoder", title:"Deep-Space Seed Stress Decoder",
          desc:"Seed stress decoding", url:"https://dr-richard-barker.github.io/deepspace-seed-stress-decoder/" },
        { id:"B_rappa_LLGCSS", title:"B. rapa — Scent × Radiation",
          desc:"Floral scent under radiation", url:"https://dr-richard-barker.github.io/B_rappa_LLGCSS/" },
      ]
    },
    {
      name: "Tropism & morphology",
      items: [
        { id:"Tropism_autodecoder_2026", title:"Tropism Autodecoder",
          desc:"Auto-decoder atlas of tropism", url:"https://dr-richard-barker.github.io/Tropism_autodecoder_2026/" },
        { id:"tropism-autodecoder-webtool", title:"Tropism Web Tool",
          desc:"Upload-and-decode browser tool", url:"https://dr-richard-barker.github.io/tropism-autodecoder-webtool/" },
        { id:"astroroot", title:"AstroRoot",
          desc:"Root architecture in microgravity", url:"https://dr-richard-barker.github.io/astroroot/" },
      ]
    },
    {
      name: "VEGGIE & crops",
      items: [
        { id:"VEGGIE_Tom_Red_Blue_Leaves_and_adv_roots", title:"VEGGIE Tomato (VEG-05)",
          desc:"Red/blue leaves & adventitious roots", url:"https://dr-richard-barker.github.io/VEGGIE_Tom_Red_Blue_Leaves_and_adv_roots/" },
        { id:"DeepLearning_VEG05", title:"PhysioSpace VEG-05",
          desc:"Sibling interaction + PhysioSpace", url:"https://dr-richard-barker.github.io/DeepLearning_VEG05/" },
        { id:"smallRNAseq-DREAM", title:"smallRNAseq-DREAM",
          desc:"microRNA-seq consolidation", url:"https://dr-richard-barker.github.io/smallRNAseq-DREAM/" },
      ]
    },
    {
      name: "Education & other",
      items: [
        { id:"Space_Biology_Education.io", title:"Space Biology Education",
          desc:"Teaching & outreach hub", url:"https://dr-richard-barker.github.io/Space_Biology_Education.io/" },
        { id:"bloodbowl", title:"Brutal Bowl",
          desc:"LLM game-theory arena game", url:"https://dr-richard-barker.github.io/Training_LLM_game-theory_using_bloodbowl/" },
      ]
    }
  ]
};
