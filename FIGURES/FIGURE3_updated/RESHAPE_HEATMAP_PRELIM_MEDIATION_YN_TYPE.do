**Mediation results with 144 proteins: reshaping of file**
**Copy paste from excel, each of FA, MD, ICVF, ISOVF, and OD results*
**Run the following reshape command and then send to Yi-Han to do the heat map**



cd "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER8E_ADPRSPROTDEM\MANUSCRIPT_FINAL\FIGURES\FIGURE3_updated"


**NO MEDIATION FINDINGS**" 
use NO_MEDIATION_GROUPA,clear

destring p, replace

save NO_MEDIATION_GROUPAfin,replace


reshape wide beta se z p lcl ucl, i(protein) j(fourwaydecomp, string)

save NO_MEDIATION_GROUPA_wide, replace


**INCONSISTENT MEDIATION FINDINGS**
use INCONSISTENT_MEDIATION_GROUPB,clear

replace p="0" if p=="<0.001"
destring p, replace

save INCONSISTENT_MEDIATION_GROUPBfin, replace

reshape wide beta se z p lcl ucl, i(protein) j(fourwaydecomp, string)

save INCONSISTENT_MEDIATION_GROUPB_wide, replace

**CONSISTENT MEDIATION FINDINGS**
use CONSISTENT_MEDIATION_GROUPC,clear

replace p="0" if p=="<0.001"
destring p, replace

save CONSISTENT_MEDIATION_GROUPCfin, replace


reshape wide beta se z p lcl ucl, i(protein) j(fourwaydecomp, string)

save CONSISTENT_MEDIATION_GROUPC_wide, replace

