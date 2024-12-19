* Cagé, Le Pennec, Mougin (2023)

global dir "C:\Users\Admin\Desktop\Replication Package\"

set matsize 800
 
********************************************************************************
* Set path to working directories:
********************************************************************************

global output "$dir/data/output"
global temp "$dir/data/temp"
global main "$dir/figures_tables/main_paper"
global appendix "$dir/figures_tables/appendix"

********************************************************************************
* Install ados / packages:
********************************************************************************

cap ssc install estout
cap ssc install esttab
cap ssc install eststo
cap ssc install estadd
cap ssc install estpost
cap ssc install cibar
cap ssc install ftools
cap ssc install reghdfe
cap ssc install reclink
cap ssc install coefplot

********************************************************************************
* Preliminaries

clear all

use "$output/data_cand_8897", clear

**define sets of controls 
global census90="pasdedip1990  sup1990 agri1990 ouvr1990  pop_65_plus1990  pop_15_241990"
global urbancontrol="cheflieudep nb_villes"
global cc="CC_charges_fonctio CC_produits_fonctio"
global dads = "DADS_nbestab DADS_sumwage DADS_share_top1  DADS_nbworker"
global elec = "limit_cst change_limit_cst inscrits1 elec_margin_c_88 secround_circo_c_88"

//time-varying controls only
global cand_cont = "female rerun Dincumbent mayor other_mandate"
global circo_cont="limit_cst c_female c_rerun c_Dincumbent c_mayor c_other_mandate inscrits1 nb_party1 nb_party2 nb_party3 nb_party4 nb_party5 nb_party6 nb_party7 nb_party8 nb_party9 nb_party10 nb_party11 nb_party12 nb_party13 nb_party14 nb_party15 nb_party16 nb_party17 nb_party18 DADS_nbestab DADS_sumwage DADS_share_top1  DADS_nbworker CC_charges_fonctio CC_produits_fonctio"

global trends = "ratio_local1_partydep_8188 ratio_local1_partydep_7881 ratio_local1_partydep_7378 ratio_local1_partydep_6873 ratio_local1_partydep_6768"

**fixed effects
egen id_circo=group(codegeo) 
egen id_indiv=group(codegeo party id_cand)	
egen id_yearparty=group(party year) 
egen id_year=group(year)
egen id_circoyear=group(codegeo year) 
egen id_partycirc=group(codegeo party) 

**predicted amount of Firm donations per voter in 1993
capture drop temp_* disp_*
foreach x of varlist *_i_93 *_c_93 change_limit_cst elec_margin_c_88 secround_circo_c_88{
gen temp_`x'=`x'
gen disp_`x'=`x'==.
replace temp_`x'=-1 if `x'==.
}

reghdfe dons_firms temp_* disp_* if year==1993, absorb(id_yearparty) cluster(id_circo)
predict pred93 if year==1993
drop temp_* disp_*

byso id_cand: egen pred_dons_firms=min(pred93) //within candidate
gen inter_pred_dons_firms=pred_dons_firms*id_year

byso id_partycirc: egen pred_dons_firms_p=mean(pred93) if party!="other" //within party*district
gen inter_pred_dons_firms_p=pred_dons_firms_p*id_year

**same, including initial contributions in 93
capture drop temp_* disp_*
foreach x of varlist dons_indiv_93 party_contrib_93 personal_contrib_93 *_i_93 *_c_93 change_limit_cst elec_margin_c_88 secround_circo_c_88{
gen temp_`x'=`x'
gen disp_`x'=`x'==.
replace temp_`x'=-1 if `x'==.
}

reghdfe dons_firms temp_* disp_* if year==1993, absorb(id_yearparty) cluster(id_circo)
predict pred93_2 if year==1993
drop temp_* disp_*

byso id_cand: egen pred_dons_firms_2=min(pred93_2) //within candidate
gen inter_pred_dons_firms_2=pred_dons_firms_2*id_year

byso id_partycirc: egen pred_dons_firms_p_2=mean(pred93_2) if party!="other" //within party*district
gen inter_pred_dons_firms_p_2=pred_dons_firms_p_2*id_year

save "$temp/analysis", replace

********************************************************************************
* Descriptive Statistics
* Table 1: campaign spendings and revenues
use "$output/data_cand_8897", clear

keep if year==1993 | year==1997

keep year codegeo candidat id_cand private_donation_indiv_cst private_donation_firms_cst total_expenditures_cst total_revenues_cst party_contribution_cst kind_contribution_cst other_cst personal_contribution_cst other_kind_cst sh_dons_firms sh_dons_indiv sh_party_contrib    sh_personal_contrib nb_elec

gen x=_n
reshape wide private_donation_indiv_cst private_donation_firms_cst total_expenditures_cst total_revenues_cst party_contribution_cst kind_contribution_cst other_cst personal_contribution_cst other_kind_cst sh_dons_firms sh_dons_indiv sh_party_contrib  sh_personal_contrib 	, i(x) j(year)

********************************************************************************
* Labeling
* Essential labels for Tables 1 and main analysis
cap label var total_expenditures_cst1993 "1993"
cap label var total_expenditures_cst1997 "1997"
cap label var total_revenues_cst1993 "1993"
cap label var total_revenues_cst1997 "1997"
cap label var sh_dons_firms1993 "1993"
cap label var sh_dons_firms1997 "1997"
cap label var sh_dons_indiv1993 "1993"
cap label var sh_dons_indiv1997 "1997"
cap label var sh_personal_contrib1993 "1993"
cap label var sh_personal_contrib1997 "1997"
cap label var sh_party_contrib1993 "1993"
cap label var sh_party_contrib1997 "1997"

* Labels for policy topics analysis
cap label var economy2_prob1 "Prevalence of economic issues"
cap label var social2_prob1 "Prevalence of social issues"
cap label var homeland2_prob1 "Prevalence of homeland and administration"
cap label var foreign2_prob1 "Prevalence of foreign issues"

* Labels for candidate characteristics
cap label var female "Female"
cap label var rerun "Re-run"
cap label var Dincumbent "Incumbent"
cap label var mayor "Mayor"
cap label var other_mandate "Other mandates"

* Labels for treatment variables
cap label var std_dons_firms "Firm donations"

* Labels for selection analysis
cap label var rerun_select "Candidate runs in next election"
cap label var party_select "Party runs in next election"
********************************************************************************
global 	stat_revenues total_expenditures_cst1993 total_expenditures_cst1997  ///
total_revenues_cst1993 total_revenues_cst1997  ///
sh_dons_firms1993 sh_dons_firms1997  ///
sh_dons_indiv1993 sh_dons_indiv1997  ///
sh_personal_contrib1993 sh_personal_contrib1997  ///
sh_party_contrib1993 sh_party_contrib1997 

eststo clear
estpost sum $stat_revenues, d

eststo : esttab using "$main/Table1.tex", label replace ///
cells("mean(fmt(%15.0fc %15.0fc %15.0fc %15.0fc %15.2fc %15.2fc %15.2fc %15.2fc) label(Mean)) p50(fmt(%15.0fc  %15.0fc %15.0fc %15.0fc %15.2fc %15.2fc %15.2fc %15.2fc) label(Median)) sd(fmt(%15.0fc %15.0fc %15.0fc %15.0fc %15.2fc %15.2fc %15.2fc %15.2fc) label(sd)) min(fmt(%15.0fc) label(Min)) max(fmt(%15.0fc) label(Max)) count(fmt(%15.0fc) label(N))") ///
style(tex) nonum noobs ///
mlabels("Spending (cst \euro)", pattern(1 0 0 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) ///
refcat(total_expenditures_cst1993 "\textbf{Total spending per candidate}" total_revenues_cst1993 "\textbf{Total revenues}" sh_dons_firms1993 "\textbf{Share firm donations}"sh_dons_indiv1993 "\textbf{Share individual donations}" sh_personal_contrib1993 "\textbf{Share personal contributions}" sh_party_contrib1993 "\textbf{Share party contributions}", nolabel)

********************************************************************************
* Figure 1: Revenues and firm donations across parties

use "$output/data_cand_8897", clear

keep if year==1993 | year==1997
replace private_donation_firms_cst=0 if private_donation_firms_cst==.

count if total_revenues_cst==.

gen nb_cand_wd=1
replace nb_cand_wd=. if private_donation_firms_cst==.

replace bigparty="Far-right" if bigparty=="FN"
replace bigparty="Right" if bigparty=="UMP"
replace bigparty="Socialist" if bigparty=="PS"
replace bigparty="Communist" if bigparty=="PC"
replace bigparty="Green" if bigparty=="Verts"
replace bigparty="Other" if bigparty=="other"

gen a=1
collapse (mean)   total_revenues_cst private_donation_firms_cst (sum)  sum_nb_cand_wd=nb_cand_wd sum_total_revenues_cst= total_revenues_cst sum_private_donation_firms_cst=private_donation_firms_cst a , by(party bigparty year)

preserve
gen sh_don=private_donation_firms_cst/total_revenues_cst
bys bigparty: sum sh_don
gen order=.

replace order=1 if bigparty=="Communist"
replace order=2 if bigparty=="Green"
replace order=3 if bigparty=="Socialist"
replace order=4 if bigparty=="Right"
replace order=5 if bigparty=="Far-right"
replace order=6 if bigparty=="Other"

graph bar total_revenues_cst private_donation_firms_cst , ///
stack over(year , label(labsize(vsmall)))  over(bigparty , sort(order) label(labsize(small))) /// 
bar(1, fcolor(black) lcolor(black)) bar(2, fcolor(white) lcolor(black)) ///
ytitle("Revenues in €")  ///
legend(label(1 "Other revenues")  label(2 "Firm donations") size(small)) ///
ylab(0 "0" 20000 "20,000" 40000 "40,000" 60000 "60,000" 80000 "80,000" 100000 "100,000" 120000 "120,000", labsize(small) angle(horizontal)) ///
graphregion(color(white))

graph export  "$main/Figure2.pdf", replace	

restore

********************************************************************************
* Table 6: Impact of firm donations on broad policy topics

estimates clear

local j=0

foreach out in economy2_prob1 social2_prob1 homeland2_prob1 foreign2_prob1{

local j=`j'+1

use "$temp/analysis", clear

keep if `out'!=. & sample_rest==1 & sample_did==1

* control
foreach x in $cand_cont{
gen disp_`x'=`x'==.
replace `x'=-1 if `x'==.
}

* fixed effects
byso id_indiv: gen temp=_N
keep if temp>1
drop temp

byso id_yearparty: gen temp=_N
replace id_yearparty=10000+year if temp<2 
drop temp

byso id_yearparty: gen temp=_N
keep if temp>=2 
drop temp

replace std_dons_firms=-std_dons_firms
label var std_dons_firms "Firm donations (loss)"

reghdfe `out' disp_* $cand_cont inter_pred_dons_firms std_dons_firms, absorb(id_indiv id_yearparty) cluster(id_circo)
sum `out' if year==1993
estadd scalar ymean=r(mean)
estimate store topic_`j'
}

esttab topic_* using "$main/Table6.tex", ///
replace keep(std_dons_firms)  ///
b(3) se lab nomtitles ///
unstack style(tex) lines compress star(* 0.10 ** 0.05 *** 0.01) nonotes ///
mgroups("\shortstack{Economic\\policy}" "\shortstack{Social\\policy}" "\shortstack{Homeland and\\administration}" "\shortstack{Foreign\\policy}", pattern(1 1 1 1)  prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) ///
scalars( "ymean Mean outcome before ban" "r2_within R2-Within") 

********************************************************************************
* Alternative Approach
* Table 6: Impact of firm donations on broad policy topics (didregress)
estimates clear
local j=0
foreach out in economy2_prob1 social2_prob1 homeland2_prob1 foreign2_prob1 {
    local j=`j'+1
    use "$temp/analysis", clear
    
    * Keep relevant sample
    keep if `out'!=. & sample_rest==1 & sample_did==1
    
    * Keep only candidates running twice
    byso id_indiv: gen temp=_N
    keep if temp>1
    drop temp
    
    * Generate treatment variable
    gen treat = std_dons_firms
    replace treat = -treat  // Multiply by -1 to match original interpretation
    
    * Generate post indicator
    gen post = (year == 1997)
    
    * Create control interactions with post
    foreach var in female rerun Dincumbent mayor other_mandate {
        gen v1_`var' = `var'
        replace v1_`var' = -1 if v1_`var' == .
        gen v1_`var'_post = v1_`var' * post
    }
    
    * Party-year interactions
    levelsof party, local(parties)
    foreach p of local parties {
        gen party_`p' = (party=="`p'")
        gen party_`p'_post = party_`p' * post
    }
    
    * DiD regression and store results directly
    didregress (`out' v1_*_post party_*_post inter_pred_dons_firms) ///
               (treat, continuous), ///
               group(id_indiv) time(year) vce(cluster id_circo)
    estimates store topic_`j'
    
    * Calculate R-squared within
    predict xb, xb
    egen mean_y = mean(`out')
    gen y_dm = `out' - mean_y
    egen mean_xb = mean(xb)
    gen xb_dm = xb - mean_xb
    quietly correlate y_dm xb_dm
    estadd scalar r2_within = r(rho)^2
    
    * Store mean outcome
    sum `out' if year==1993
    estadd scalar ymean = r(mean)
    
    * Clean up
    drop xb mean_y y_dm mean_xb xb_dm
}

* Output table
esttab topic_* using "$main/Table6_didregress.tex", ///
    replace keep(treat) ///
    coeflabel(treat "Firm donations (loss)") ///
    b(3) se lab nomtitles ///
    unstack style(tex) lines compress star(* 0.10 ** 0.05 *** 0.01) nonotes ///
    mgroups("\shortstack{Economic\\policy}" "\shortstack{Social\\policy}" "\shortstack{Homeland and\\administration}" "\shortstack{Foreign\\policy}", ///
    pattern(1 1 1 1) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) ///
    scalars("ymean Mean outcome before ban" "r2_within R2-Within")

	
********************************************************************************	
* Post-Treatment Heterogeneity Test - Subgroup Analysis
estimates clear
local j = 0

foreach het in female Dincumbent mayor {
    * For each value of the heterogeneity variable (0/1)
    forvalues val = 0/1 {
        local j = `j' + 1
        use "$temp/analysis", clear
        
        * Keep relevant sample and subgroup
        keep if economy2_prob1!=. & sample_rest==1 & sample_did==1
        keep if `het' == `val'
        
        * Keep only candidates running twice
        byso id_indiv: gen temp=_N
        keep if temp>1
        drop temp
        
        * Handle missing values in controls
        foreach x in $cand_cont {
            gen disp_`x'=`x'==.
            replace `x'=-1 if `x'==.
        }
        
        * Prepare party-year groups
        byso id_yearparty: gen temp=_N
        replace id_yearparty=10000+year if temp<2 
        drop temp
        
        byso id_yearparty: gen temp=_N
        keep if temp>=2 
        drop temp
        
        * Flip sign of treatment for interpretation
        replace std_dons_firms=-std_dons_firms
        label var std_dons_firms "Firm donations (loss)"
        
        * Run reghdfe for subgroup
        reghdfe economy2_prob1 disp_* $cand_cont inter_pred_dons_firms std_dons_firms, ///
            absorb(id_indiv id_yearparty) cluster(id_circo)
            
        * Store results
        sum economy2_prob1 if year==1993
        estadd scalar ymean=r(mean)
        estimate store hetero_`het'_`val'
    }
}

* Output table
esttab hetero_* using "$main/Table6_heterogeneity.tex", ///
    replace keep(std_dons_firms) ///
    b(3) se lab nomtitles ///
    unstack style(tex) lines compress star(* 0.10 ** 0.05 *** 0.01) nonotes ///
    mgroups("\shortstack{Male\\Female}" "\shortstack{Non-Incumbent\\Incumbent}" "\shortstack{Non-Mayor\\Mayor}", ///
    pattern(1 1) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) ///
    scalars("ymean Mean outcome before ban" "r2_within R2-Within")
	
********************************************************************************	
* Appendix: Robustness checks
* Table E.11: Firm donations and selection into sample
estimates clear 

use "$temp/analysis", clear

keep if sample_rest==1 & sample_did==1

gen temp=text_size_sr1!=.
byso id_indiv: egen temp2=total(temp)
gen nomissing_manif=temp2==2
drop temp*
ta nomissing_manif if year==1993

gen temp=montant_total!=.
byso id_indiv: egen temp2=total(temp)
gen nomissing_dons=temp2==2
drop temp*
ta nomissing_dons if year==1993

* control
foreach x in $cand_cont{
gen disp_`x'=`x'==.
replace `x'=-1 if `x'==.
}

byso id_yearparty: gen temp=_N
replace id_yearparty=10000+year if temp<2 
drop temp

byso id_yearparty: gen temp=_N
keep if temp>=2 
drop temp

keep if year==1993

reghdfe rerun_select  disp_* $cand_cont pred_dons_firms std_dons_firms, absorb(id_yearparty) cluster(id_circo)	
sum rerun_select if year==1993
estadd scalar ymean=r(mean)
estimate store select_1

reghdfe nomissing_manif disp_* $cand_cont pred_dons_firms std_dons_firms if rerun_select==1, absorb(id_yearparty) cluster(id_circo)	
sum nomissing_manif if year==1993
estadd scalar ymean=r(mean)
estimate store select_2


use "$temp/analysis", clear

keep if sample_rest==1 & sample_did==1

keep if bigparty!="other" & bigparty!=""

foreach x in $cand_cont{
byso id_partycirc year: egen p_`x'=mean(`x')
gen disp_`x'=p_`x'==.
replace p_`x'=-1 if p_`x'==.
}

byso id_yearparty: gen temp=_N
replace id_yearparty=10000+year if temp<2 
drop temp

byso id_yearparty: gen temp=_N
keep if temp>=2 
drop temp

keep if year==1993

reghdfe party_select  disp_* $cand_cont pred_dons_firms std_dons_firms, absorb(id_yearparty) cluster(id_circo)
sum party_select if year==1993
estadd scalar ymean=r(mean)
estimate store select_3

esttab select_*  using "$appendix/TableE11.tex", ///
replace keep(std_dons_firms) b(3) se lab nomtitles ///
style(tex) lines compress star(* 0.10 ** 0.05 *** 0.01) nonotes ///
mgroups("\shortstack{Candidate in\\next election}" "\shortstack{Manifesto\\available}" "\shortstack{Party in\\next election}", pattern(1 1 1)  prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) ///
scalars("ymean Mean outcome before ban" "r2_within R2-Within")	