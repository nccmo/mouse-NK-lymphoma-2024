mNK_SG <- c("mNK_4285_SG", "mNK_4284_SG", "mNK_3714_SG", "mNK_3856_SG", "mNK_4214_SG")
mNK_SP <- c("mNK_2754_SP", "mNK_3660_SP", "mNK_4285_SP", "mNK_4284_SP", 
            "mNK_3714_SP", "mNK_3856_SP", "mNK_4453_SP", "mNK_4457_SP", "mNK_3777_SP", "mNK_3797_SP", "mNK_4214_SP")
mNK_Tu <- c("mNK_4453_Tu", "mNK_4457_Tu")

disease_type_mNK_level = c("mNK_SG", "mNK_SP", "mNK_Tu")

disease_type2_mNK_level = c("CreN_pre", "CreP_pre", "CreP_Tumor", "CreLMP1P_Tumor")

disease_type3_NKmouse_level = c("SG_CreN_pre", "SP_CreN_pre", 
                                "SG_CreP_pre", "SP_CreP_pre", 
                                "SG_CreP_Tumor", "SP_CreP_Tumor", 
                                "Tu_CreLMP1P_Tumor", "SP_CreLMP1P_Tumor")

disease_type4_NKmouse_level = c("SG_pre", "SP_pre", 
                            "SG_CreP_Tumor", "SP_CreP_Tumor", 
                            "SG_CreLMP1P_Tumor", "CreLMP1P_Tumor")

disease_type5_NKmouse_level = c("SG_pre", "SP_pre", 
                                "SG_CreP_Tumor", "SP_CreP_Tumor", 
                                "SG_CreLMP1P_Tumor", "SP_CreLMP1P_Tumor", "Tu_CreLMP1P_Tumor")

disease_type6_NKmouse_level = c("mNK_SG_Tu", "mNK_SP")

disease_type7_NKmouse_level = c("Cre_pre", "CreP_Tumor", "CreLMP1P_Tumor")

annotationSeurat1_NKmouse_filt2_level = c("NM01_T_NK01", "NM02_T_NK02", "NM03_T_NK03", 
                                          "NM04_B01", "NM05_B02", "NM06_B03", "NM07_B04", 
                                          "NM08_Myeloid01",
                                          "NM09_Myeloid02", "NM10_Myeloid03", "NM11_Myeloid04", 
                                          "NM12_Myeloid05",
                                          "M01_Tumor01", "M02_Tumor02", "M03_Tumor03",
                                          "M04_Tumor04", "M05_Tumor05", "M06_Tumor06",
                                          "M07_Tumor07")

annotationSeurat1_NKmouse_filt2_level_HK <- c(
  "NM01_T_NK01",
  "NM02_T_NK02",
  "NM03_T_NK03",
  "NM04_T_NK04",
  "NM05_B01",
  "NM06_B02",
  "NM07_B03",
  "NM08_B04",
  "NM09_Myeloid01",
  "NM10_Myeloid02",
  "NM11_Myeloid03",
  "NM12_Myeloid04",
  "NM13_Myeloid05",
  "M01_Tumor01",
  "M02_Tumor02",
  "M03_Tumor03",
  "M04_Tumor04",
  "M05_Tumor05",
  "M06_Tumor06",
  "M07_Tumor07"
)

celltype1_NKmouse_filt2_level = c("T_NK", 
                                  "B_immature", "B_mature", "B_Hsp70", "PBL", 
                                  "Mono", "MDSC", "Mac", "cDC", "pDC", 
                                  "Tumor_LMPn", "Tumor_LMPp")

celltype2_NKmouse_filt2_level = c("T_NK", 
                                  "B", 
                                  "Myeloid", 
                                  "Tumor")

celltype3_NKmouse_filt2_level = c("CD4T", "CD8T", "CD4CD8DN", "NK",
                                  "B", 
                                  "Myeloid", 
                                  "Tumor")

celltype4_NKmouse_filt2_level = c("CD4T", "CD8T", "CD4CD8DN", "NK", 
                                  "B", 
                                  "Myeloid", 
                                  "Tumor_LMPn", "Tumor_LMPp")

celltype4_NKmouse_filt2_level_HK = c("CD4T", "CD8T", "CD4CD8DN", 
                                  "B", 
                                  "Myeloid", 
                                  "NK",
                                  "Tumor_LMPn", "Tumor_LMPp")

celltype3_NKmouse_filt2_level_HK = c("CD4T", "CD8T", "CD4CD8DP", "CD4CD8DN", "NK",
                                     "B", 
                                     "Myeloid", 
                                     "Tumor")

annotationSeurat1_T_NK5_NKmouse_filt2_level = c("CD4T", "CD8T", "CD4CD8DN", "NK",
                                                "B_immature", "B_mature", "B_Hsp70", "PBL", 
                                                "Mono", "MDSC", "Mac", "cDC", "pDC", 
                                                "Tumor_LMPn", "Tumor_LMPp")

annotationSeurat1_T_NK5_NKmouse_filt2_level_HK = c("CD4T", "CD8T", "CD4CD8DP", "CD4CD8DN", "NK",
                                                "B_immature", "B_mature", "B_Hsp70", "PBL", 
                                                "Mono", "MDSC", "Mac", "cDC", "pDC", 
                                                "Tumor_LMPn", "Tumor_LMPp")

annotationSeurat2_NKmouse_filt2_level = c("CD4Tsub", "TFH", 
                                          "TREG", "CD8T_N",
                                          "CD8T_GZMK", 
                                          "CD4CD8DN", "NK", "ILC3",
                                          "B_immature", "B_mature", "B_Hsp70", "PBL", 
                                          "Mono", "MDSC", "Mac", "cDC", "pDC", 
                                          "Tumor_LMPn", "Tumor_LMPp")

annotationSeurat2_NKmouse_filt2_level_HK = c("CD4T_N", "CD4T_CM", "TREG", 
                                             "CD4T_other1", "CD4T_other2", 
                                             "CD8T_N", "CD8T_GZMK", 
                                             "CD4CD8DN", "NK", "ILC3",
                                             "B_immature", "B_mature", "PBL", 
                                             "Mono", "MDSC", "Mac", "cDC", "pDC", 
                                             "Tumor_LMPn", "Tumor_LMPp")

annotationSeurat2_NKmouse_filt2_level_HK2 = c("CD4T_N", "TREG", 
                                             "CD4T_other1", "CD4T_other2", 
                                             "CD8T_N", "CD8T_GZMK", 
                                             "CD4CD8DN", "ILC3",
                                             "B_immature", "B_mature", "PBL", 
                                             "Mono", "MDSC", "Mac", "cDC", "pDC", 
                                             "NK", 
                                             "Tumor_LMPn", "Tumor_LMPp")

annotationSeurat2_micetype_NKmouse_filt2_level = c("CD4Tsub_Cre_pre", "CD4Tsub_CreP_Tumor", "CD4Tsub_CreLMP1P_Tumor", 
                                                   "TFH_Cre_pre", "TFH_CreP_Tumor", "TFH_CreLMP1P_Tumor", 
                                                   "TREG_Cre_pre", "TREG_CreP_Tumor", "TREG_CreLMP1P_Tumor", 
                                                   "CD8T_N_Cre_pre", "CD8T_N_CreP_Tumor", "CD8T_N_CreLMP1P_Tumor",
                                                   "CD8T_GZMK_Cre_pre", "CD8T_GZMK_CreP_Tumor", "CD8T_GZMK_CreLMP1P_Tumor", 
                                                   "CD4CD8DN_Cre_pre", "CD4CD8DN_CreP_Tumor", "CD4CD8DN_CreLMP1P_Tumor", 
                                                   "NK_Cre_pre", "NK_CreP_Tumor", "NK_CreLMP1P_Tumor", 
                                                   "ILC3_Cre_pre", "ILC3_CreP_Tumor", "ILC3_CreLMP1P_Tumor",
                                                   "B_immature_Cre_pre", "B_immature_CreP_Tumor", "B_immature_CreLMP1P_Tumor", 
                                                   "B_mature_Cre_pre", "B_mature_CreP_Tumor", "B_mature_CreLMP1P_Tumor", 
                                                   "B_Hsp70_Cre_pre", "B_Hsp70_CreP_Tumor", "B_Hsp70_CreLMP1P_Tumor", 
                                                   "PBL_Cre_pre", "PBL_CreP_Tumor", "PBL_CreLMP1P_Tumor", 
                                                   "Mono_Cre_pre", "Mono_CreP_Tumor", "Mono_CreLMP1P_Tumor", 
                                                   "MDSC_Cre_pre", "MDSC_CreP_Tumor", "MDSC_CreLMP1P_Tumor", 
                                                   "Mac_Cre_pre", "Mac_CreP_Tumor", "Mac_CreLMP1P_Tumor", 
                                                   "cDC_Cre_pre", "cDC_CreP_Tumor", "cDC_CreLMP1P_Tumor", 
                                                   "pDC_Cre_pre", "pDC_CreP_Tumor", "pDC_CreLMP1P_Tumor", 
                                                   "Tumor_LMPn_Cre_pre", "Tumor_LMPn_CreP_Tumor", "Tumor_LMPn_CreLMP1P_Tumor", 
                                                   "Tumor_LMPp_Cre_pre", "Tumor_LMPp_CreP_Tumor", "Tumor_LMPp_CreLMP1P_Tumor")

annotationSeurat2_micetype2_NKmouse_filt2_level = c("CD4Tsub_pre", "CD4Tsub_Tumor", 
                                                   "TFH_pre", "TFH_Tumor", 
                                                   "TREG_pre", "TREG_Tumor", 
                                                   "CD8T_N_pre", "CD8T_N_Tumor",
                                                   "CD8T_GZMK_pre", "CD8T_GZMK_Tumor", 
                                                   "CD4CD8DN_pre", "CD4CD8DN_Tumor", 
                                                   "NK_pre", "NK_Tumor", 
                                                   "ILC3_pre", "ILC3_Tumor",
                                                   "B_immature_pre", "B_immature_Tumor", 
                                                   "B_mature_pre", "B_mature_Tumor", 
                                                   "B_Hsp70_pre", "B_Hsp70_Tumor", 
                                                   "PBL_pre", "PBL_Tumor", 
                                                   "Mono_pre", "Mono_Tumor", 
                                                   "MDSC_pre", "MDSC_Tumor", 
                                                   "Mac_pre", "Mac_Tumor", 
                                                   "cDC_pre", "cDC_Tumor", 
                                                   "pDC_pre", "pDC_Tumor", 
                                                   "Tumor_LMPn_pre", "Tumor_LMPn_Tumor", 
                                                   "Tumor_LMPp_pre", "Tumor_LMPp_Tumor")

annotationSeurat1_T_NK5v2_NKmouse_filt2_level = c("CD4T", "CD8T", "CD4CD8DN", "NK",
                                                  "B_immature", "B_mature", "B_Hsp70", "PBL", 
                                                  "Mono", "MDSC", "Mac", "cDC", "pDC", 
                                                  "M01_Tumor01", "M02_Tumor02", "M03_Tumor03",
                                                  "M04_Tumor04", "M05_Tumor05", "M06_Tumor06",
                                                  "M07_Tumor07")

celltype4_level3 = c("CD4T", "CD8T", "CD4CD8DN", "NK", "B", "Myeloid", "Tumor_LMPn", "Tumor_LMPp")
