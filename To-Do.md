# To do
- [ ] For dataset CAL (and others): investigate whether the survival event for dmfs_status is "distant metastasis" or "distant metastasis OR death". In the case that a patient is deceased and did not experience distant metastasis, what is the expected event status?
- [ ] TCGA days_to_death: for patients with status "living", populate days to last follow-up instead of NA
- [ ] When merging data sets with rbind: ensure there are factors do not get coerced into integers
