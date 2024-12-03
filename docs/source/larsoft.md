# HyperonProduction Reference

:::{card}
:link: https://github.com/playonverbs/HyperonProduction

This page is a standalone reference for <github:playonverbs/HyperonProduction>,\
a module for the LArSoft framework written in C++.
:::

::::::{card} `HyperonProduction_module.cc`

::::{cpp:struct} hyperon::Config

Contains fields that correspond to {abbr}`FHiCL (Fermilab Hierarchical
Configuration Language)` configuration parameters.

See
[this](https://indico.fnal.gov/event/11857/sessions/1051/attachments/6785/8812/LArSoftUsability_workshop_June2016_knoepfel.pdf)
for an informal FHiCL reference.

:::{cpp:member} Atom<std::string> fPandoraRecoLabel
Label for pandoraPatRec reconstruction products.
:::
:::{cpp:member} Atom<std::string> fFlashMatchRecoLabel
Label for pandora reconstruction products.
:::
:::{cpp:member} Atom<std::string> fTrackLabel
Label for `recob::Track`.
:::
:::{cpp:member} Atom<std::string> fShowerLabel
Label for `recob::Shower`
:::
:::{cpp:member} Atom<std::string> fCaloLabel
Label for `anab::Calorimetry`
:::
:::{cpp:member} Atom<std::string> fGeneratorLabel
Label for `simb::MCTruth`
:::
:::{cpp:member} Atom<std::string> fG4Label
Label for `simb::MCParticle`
:::
:::{cpp:member} Atom<std::string> fPIDLabel
Label for `anab::ParticleID`
:::
:::{cpp:member} Atom<std::string> fHitLabel
Label for `recob::Hit`
:::
:::{cpp:member} Atom<std::string> fTrackHitAssnsLabel
Label for associations between `recob::Track` and `recob::Hit`
:::
:::{cpp:member} Atom<std::string> fHitTruthAssnsLabel
Label for associations between `simb::MCParticle`, `recob::Hit` and
`anab::BackTrackerHitMatchingData`.
:::
:::{cpp:member} Atom<std::string> fPOTSummaryLabel
Label for POT Summary data.
:::
:::{cpp:member} Atom<bool> fIsData
Flag to indicate if the input is Data.
:::
:::{cpp:member} Atom<bool> fDebug = false
Flag to enable debug messages.
:::
::::

:::::{cpp:class} hyperon::HyperonProduction : art::EDAnalyzer
:::{cpp:function} void analyze(art::Event const& evt) override
:::
:::{cpp:function} void beginJob() override
:::
:::{cpp:function} void endJob() override
:::
:::{cpp:function} void beginSubRun(const art::SubRun& sr)
:::
:::{cpp:function} void endSubRun(const art::SubRun& sr)
:::
:::::
::::::

:::::{card} `util.h`

::::{cpp:enum-struct} hyperon::util::GenEventType : int
Represents the neutrino interaction type as recorded by the event generator.

:::{cpp:enumerator} QEL = 0
:::
:::{cpp:enumerator} RES = 1
:::
:::{cpp:enumerator} DIS = 2
:::
:::{cpp:enumerator} COH = 3
:::
:::{cpp:enumerator} ElectronScatter = 5
:::
:::{cpp:enumerator} MEC = 10
:::
:::{cpp:enumerator} Diffractive = 11
:::
:::{cpp:enumerator} HYP = 1095
:::
:::{cpp:enumerator} Other
:::
::::

:::{cpp:function} std::string hyperon::util::GetCCNC(int ccnc)
Returns the string representation of the charged/neutral current interaction
output by the generator.

:retval "CC": For a charged-current interaction
:retval "NC": For a neutral-current interaction
:retval "None": For neither.
:::

:::{cpp:function} std::string hyperon::util::GetEventType(int mode)
Returns the string representation of the event interaction type given by the
generator.
:::

---

:::{cpp:function} template <typename T> std::vector<art::Ptr<T>> hyperon::util::GetProductVector(const art::Event &e, const std::string &label)

:tparam T: A type that can be stored within an artroot event file.

Returns all data products of type `T` associated with the event that have a
module label `label`.
:::

:::{cpp:function} template <typename T, typename U> std::vector<art::Ptr<T>> hyperon::util::GetAssocProductVector(const art::Ptr<U> &pProd, const art::Event &e, const std::string &label, const std::string &assocLabel)
:::

:::{cpp:function} template <typename T, typename U> art::Ptr<T> hyperon::util::GetAssocProduct(const art::Ptr<U> &pProd, const art::Event &e, const std::string &label, const std::string &assocLabel)
:::

:::{cpp:function} art::Ptr<simb::MCParticle> hyperon::util::getAssocMCParticle(art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &hittruth, const std::vector<art::Ptr<recob::Hit>> &hits, float &purity)

This function back-tracks hits to simulated particles when given an
`art::FindManyP`, a list of `recob::Hit` and an associated
`ana::BackTrackerHitMatchingData`.

This function checks for which `simb::MCParticle` deposited the highest energy
that corresponds with a particular `recob::Hit`.
:::

:::{cpp:function} bool hyperon::util::posMatch(TVector3 p1, TVector3 p2, const double _epsilon = 0.0001)

Checks if two vectors are close. `_epsilon` catches floating-point comparison
issues.
:::

---

:::{cpp:var} constexpr int hyperon::pdg::Lambda = 3122
:::
:::{cpp:var} constexpr int hyperon::pdg::NeutralKaon = 311
:::
:::{cpp:var} constexpr int hyperon::pdg::SignaZero = 3212
:::

:::{cpp:function} inline bool hyperon::pdg::isHyperon(const int pdgCode)
:::
:::{cpp:function} inline bool hyperon::pdg::isHyperon(const art::Ptr<simb::MCParticle> p)
:::

:::{cpp:function} inline bool hyperon::pdg::isPion(const int pdgCode)
:::
:::{cpp:function} inline bool hyperon::pdg::isPion(const art::Ptr<simb::MCParticle> p)
:::

:::{cpp:function} inline bool hyperon::pdg::isNucleon(const int pdgCode)
:::
:::{cpp:function} inline bool hyperon::pdg::isNucleon(const art::Ptr<simb::MCParticle> p)
:::

:::{cpp:function} inline bool hyperon::pdg::isLepton(const int pdgCode)
:::
:::{cpp:function} inline bool hyperon::pdg::isLepton(const art::Ptr<simb::MCParticle> p)
:::

:::{cpp:function} inline bool hyperon::pdg::isNeutrino(const int pdgCode)
:::
:::{cpp:function} inline bool hyperon::pdg::isNeutrino(const art::Ptr<simb::MCParticle> p)
:::

:::{cpp:function} inline bool hyperon::pdg::isKaon(const int pdgCode)
:::
:::{cpp:function} inline bool hyperon::pdg::isKaon(const art::Ptr<simb::MCParticle> p)
:::

:::::
