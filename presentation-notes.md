# Speaker Notes — "Reservations and the IGM of Native Americans, 1900–1940"
*~30 minutes. Appendix slides available but not in the main flow.*

---

## Slide 1 — Title (~1 min)
- Welcome, brief orientation: this is early-stage work on a group that is almost entirely absent from the intergenerational mobility literature.
- Period: 1900–1940, the assimilation era — a particularly formative and disrupted stretch of Native American economic history.

---

## Slide 2 — Motivation and Research Question (~3 min)
- Open with the Meriam Report quote — it captures the prevailing federal view of Native American communities as self-reinforcing in their disadvantage. That framing motivates the question: was that actually true?
- The Meriam Report (1928) documented low incomes, poor health, substandard education. That was not a surprise finding — it was basically confirming what was already known.
- Fast forward: Native Americans remain one of the most economically disadvantaged demographic groups today (Chetty et al. 2018, 2020 — among the lowest-income groups in absolute upward mobility estimates).
- Two questions driving the paper:
  1. What was the *structure* of occupational persistence in the early 20th century? Not just "was there mobility" but what kind.
  2. How does the *measurement* of occupational position shape what we see? This is a methodological contribution — the coding choice matters a lot, and I'll show why.

---

## Slide 3 — Preview of Results (~1.5 min)
- High headline mobility, but it was mostly *structural* — driven by aggregate shifts in the occupational distribution, not by sons swapping places with each other.
- The dominant story: sons of farmers became manual wage laborers. That's not upward mobility — it's displacement.
- The top of the distribution (nonmanual) was rare and didn't reproduce itself.
- Nonemployment was sticky — sons of nonemployed fathers were much more likely to be nonemployed themselves.
- Enormous regional variation: the South barely moved; the North proletarianized almost completely; Oklahoma is its own story.

---

## Slide 4 — Historical and Economic Context (~2.5 min)
- The Dawes Act (1887) is the key structural shock. Allotment broke up communal landholdings — tribal land fell from roughly 138 million to 48 million acres by 1934. That's a 65% reduction. Carlson (1981) and Dippel et al. (2024) document that allotment actually *reduced* agricultural productivity.
- The Meriam Report's finding that reservation land was unsuitable for family farming matters here — what was left after allotment often wasn't farmable, which directly shapes the economic transitions we see.
- Federal management misaligned incentives for labor force participation — this shows up in the nonemployment numbers.
- Indian Reorganization Act (1934) ended allotment and is essentially the close of this study's window. We're capturing the tail end of a major disruption.
- Don't spend too long here — the point is just that this wasn't a stable period. The economy was being actively restructured from above.

---

## Slide 5 — Data (~2.5 min)
- Three sources: IPUMS full-count census data 1900–1940, Census Linking Project crosswalks (Tian et al. 2023), and Pritzker (2000) for tribal/regional classification.
- Sample: Native American men ages 20–44 in 1940, linked to fathers in earlier censuses.
- **N = 12,246 father-son pairs** — that's 22.6% of the 1940 AIAN male census population that got linked. Linkage rates are a known issue with census-linking approaches, and I handle it via propensity score weighting.
- Weighting target: ATC — I want the linked sample to look like the *full* 1940 AIAN male population, not just the linked subsample. After trimming at the 99th weight percentile, the effective N is ~8,300 (68% of the linked sample). Balance is good: post-trimming max |SMD| is 0.025.
- If asked: the appendix has the full weighting details and the linkage scenarios.

---

## Slide 6 — Occupational Groupings (~1.5 min)
- Two-level hierarchy: **macro** (4 categories) and **meso** (6 categories).
- Macro: Farming, Manual, Non-manual, Not Employed.
- Meso breaks manual into crafts/gov't services vs. semi-skilled/unskilled, and farming into farmers vs. farmworkers.
- The macro level is the main unit of analysis. Meso is in the appendix for completeness.
- Brief note on why this grouping: it follows established practice (Durlauf et al. 2024 use a similar scheme for white and Black men, which I use as a benchmark), and it's appropriate given the occupational distribution of this population — there aren't enough nonmanual workers to subdivide that category meaningfully.

---

## Slide 7 — Occupation Classification: Main vs. Attachment-Adjusted (~2 min)
- This is a methodological contribution. The question: how do you handle fathers observed across multiple census years, and what do you do with men who are in the labor force but currently unemployed?
- **Main spec**: Take the modal occupation across census years, preferring any employed code if one exists. Son's 1940 code at face value.
- **Attachment-adjusted**: No employment preference for modal; if the man was not in the labor force or unemployed, classify him as nonemployed — even if he has a reported occupation code.
- This reclassifies **10% of fathers and 15% of sons**.
- The consequence: nonemployment rises from 6% to 15% for fathers, and from 12% to **26%** for sons. That's not a small difference — the face-value codes significantly understate labor market detachment.
- I'll show both specs throughout. The attachment-adjusted is arguably more accurate for this population.

---

## Slide 8 — Descriptive Data (~1.5 min)
- Fathers (main spec): 66% farming, 23% manual, 6% nonemployed, 5% nonmanual. Average age ~62.
- Sons (main spec): 40% farming, 43% manual, 12% nonemployed, 5% nonmanual. Average age ~30.
- That shift — from 66% farming fathers to 43% manual sons, with farming dropping to 40% — is the dominant aggregate trend. The occupational structure changed dramatically in one generation.
- Three regions make up 60% of the sample: Southwest (24%), Oklahoma (20%), Plains (17%).

---

## Slide 9 — Markov Model of Mobility (~1.5 min)
- Quick methods note — keep this brief, the audience will either know this or doesn't need the math.
- Transition matrix: rows are father's occupation, columns are son's occupation. Each row sums to 1 — so each entry is a conditional probability.
- Steady state: what does the occupational distribution converge to if we apply the matrix indefinitely? It's a diagnostic for long-run equilibrium, not a prediction.
- The mobility measures (next slide) all derive from this matrix and the observed distribution.

---

## Slide 10 — Transition Probabilities, Macro (~3 min)
- Walk the audience through the key cells. Don't read every number — pick the story.
- **Farming row**: Most sons of farmers stayed in farming or entered manual. Very few went nonmanual. The farming-to-farming diagonal is high, but the farming-to-manual cell is also large.
- **Nonmanual row**: Sons of nonmanual fathers had a reasonably high chance of staying nonmanual, but the class was too small to matter much in aggregate.
- **Nonemployed row**: Sons of nonemployed fathers had elevated nonemployment rates themselves — persistence at the bottom.
- **Big takeaway**: the matrix is heavily downward-triangular. Upward mobility to nonmanual was rare regardless of starting point.
- Briefly note the appendix button for White and Black comparisons (Durlauf et al. 2024) — benchmarking will come later in revisions.

---

## Slide 11 — Transition Probabilities, Attachment-Adjusted (~1 min)
- The nonemp diagonal is now much more visible. When you don't force men with weak labor force attachment into employed occupation categories, the nonemployment persistence becomes stark.
- The farming-to-farming probability shrinks — some fathers who looked like farmers under the main spec are now classified as nonemployed, redistributing that mass.
- The two specs tell similar stories about direction but different stories about magnitude of nonemployment entrenchment.

---

## Slide 12 — Mobility Definitions (~1 min)
- Quick: don't linger.
- **Overall mobility (OM)**: 1 minus the weighted probability of staying in the same category. Simple.
- **Structural mobility (SM)**: The minimum mobility required given the change in the aggregate distribution. This is "forced" mobility — the economy changed, so people had to move even if there was no rank-switching.
- **Exchange mobility (EM)**: What's left over. True rank-switching — sons and fathers swapping places in a stable distribution. EM = OM − SM.
- The decomposition matters because high OM can look impressive but if it's all SM, it just means the distribution shifted, not that anyone got ahead on merit or connections.

---

## Slide 13 — OM/EM/SM Figure (~2 min)
- Nationally: overall mobility is around 0.50 — half of sons ended up in a different category than their fathers.
- But the SM share dominates EM. Most of that mobility was structural — forced by the collapse of the farming sector and the expansion of manual wage labor.
- Exchange mobility was lower, meaning relatively little rank-switching. The relative standings of families in the distribution were more persistent than the headline number suggests.
- Attachment-adjusted: both OM and the breakdown shift, but the structural-dominant story holds.
- **Punchline**: high gross mobility, low exchange mobility — a pattern consistent with a population being economically displaced, not one climbing the ladder.

---

## Slide 14 — Section Break: Measures by Region (~15 sec)
- Brief: regional heterogeneity is one of the main empirical findings. Five maps follow.

---

## Slide 15 — Structural Mobility by Region (~1.5 min)
- SM varies substantially. Regions that saw the largest occupational shifts — particularly toward manual labor — show high SM.
- Great Lakes and Northwest are high-SM regions: the occupational distribution shifted dramatically in one generation. Sons of farmers almost uniformly entered manual wage labor.
- Southwest is lower — there, farming persisted more.
- The map makes the geographic clustering visible: the Northern regions were more economically disrupted.

---

## Slide 16 — Exchange Mobility by Region (~1.5 min)
- EM is universally low relative to SM, but varies.
- Oklahoma stands out — it has higher EM than most regions. Consistent with the narrative of Oklahoma as economically dynamic (oil, commerce), where there were more actual opportunities for rank-switching.
- South has low EM: the distribution barely changed *and* there was little rank-switching. Near-total stasis.
- Regional numbers: South has OM = 0.37, P(up) = 0.06, P(down) = 0.44. Northwest: OM = 0.67, P(down) = 0.81. Great Lakes: OM = 0.64, P(down) = 0.83.

---

## Slide 17 — P(Farming | Father Farming) by Region (~1.5 min)
- How sticky was farming? Very sticky in the South and Southwest. If your father was a farmer, you were very likely to remain a farmer.
- In the Northern regions (Great Lakes, Northwest, Northeast), the farming-to-farming probability collapses. Sons of farmers there almost uniformly exited farming — mostly into manual wage labor.
- This is the proletarianization story. It wasn't voluntary upward mobility; it was the dissolution of the farm economy.

---

## Slide 18 — P(Manual | Father Farming) by Region (~1.5 min)
- The flip side of the previous map. Where farming persistence is low, the farming-to-manual probability is high.
- Great Lakes and Northwest: very high probability that a son of a farmer became a manual laborer. Near-complete proletarianization in those regions.
- This is the dominant mobility pathway nationally — it accounts for most of the overall mobility in the transition matrix.

---

## Slide 19 — P(Nonemp | Father Nonemp) by Region (~1.5 min)
- Nonemployment persistence is high nearly everywhere, but especially in Oklahoma and some Plains regions.
- Oklahoma's combination — moderate EM but high nonemployment persistence — is interesting. There were opportunities for those who were employed, but the nonemployed were stuck.
- Attachment-adjusted (appendix): nonemp persistence goes up everywhere. Once you account for weak labor force attachment, the bottom of the distribution is more entrenched than the face-value codes suggest.

---

## Slide 20 — Conclusion (~3 min)
- Three takeaways:
  1. **High overall mobility was mostly structural.** The farming sector collapsed; sons didn't climb the ladder, they were displaced from farming into manual wage labor. This is important for interpreting optimistic mobility statistics about this period.
  2. **Nonemployment is understated by face-value codes.** The attachment-adjusted spec reveals a more entrenched bottom of the distribution. Methodologically, this matters for any study using historical occupation codes for populations with high labor market detachment.
  3. **Regional heterogeneity is large and not yet explained.** The South is different from the North is different from Oklahoma. The next step is to explain *why* — reservation presence, land quality, proximity to labor markets, etc.

- Open questions (be ready for these):
  - Which spec is "right"? The attachment-adjusted is arguably more accurate but requires a judgment call about what "nonemployment" means in this context.
  - Is nonemployment persistence a labor market outcome (discrimination, skills) or a structural feature of the reservation economy (no jobs to take)? The data can't fully distinguish these — that's a limitation.

---

*Appendix slides (weighting details, linkage scenarios, meso matrices, alt-spec maps) available for Q&A.*
