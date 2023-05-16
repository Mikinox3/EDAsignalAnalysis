import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
import scikit_posthocs as sp
import numpy as np

'''
biblio :
from https://www.reneshbedre.com/blog/anova.html
https://pythonawesome.com/multiple-pairwise-comparisons-post-hoc-tests-in-python/

'''
def Anova(NTval,TBval,TSAval,df2):
    # stats f_oneway functions takes the groups as input and returns ANOVA F and p value
    fvalue, pvalue = stats.f_oneway(NTval.dropna(),TBval.dropna(), TSAval.dropna())
    print("ANOVA : ",fvalue, pvalue)

    # Ordinary Least Squares (OLS) model
    model = ols('Value ~ C(Groupe)', data=df2).fit()
    anova_table = sm.stats.anova_lm(model, typ=2) 
    print(anova_table)
    return pvalue


#Stage 0 set the parameters
#############################
Plot= True
 
#Stage 1 : get the file
##############################
aFile = 'C:/Users/mathi/Documents/Visual code/M1 Projet Rembau/MRiseTime2.csv'
df = pd.read_csv(aFile,sep=';', header=0)

#Stage 2 : select what you want to compare between groups in ANOVA 1W or KW test
##############################
#select = ['Ineg']
select = ['Ipos','Ineu','Ineg'] # group comparison for img all valences
#select = ['Fpos','Fneu','Fneg'] # group comparison for all film valences 
df2 = df[df['Valence'].isin(select)]
print(df2)

#Stage 2.1 : Plot the data to be compared
##############################
if Plot :
    ax = sns.boxplot(x='Groupe', y='Value', data=df2, color='#99c2a2')
    ax = sns.swarmplot(x='Groupe', y="Value", data=df2, color='#7d0013')
    plt.title('RMSC for +/=/- images, group comparison')
    plt.show()

#prepare the data for stats
NTval=df2.loc[df2['Groupe']=='NT',['Value']]
TBval=df2.loc[df2['Groupe']=='TB',['Value']]
TSAval=df2.loc[df2['Groupe']=='TSA',['Value']]

#Stage 3 : check the parametric test conditions 
##############################

# 3.1 independent quantitative values : CHECKED

# 3.2 normality (if n<30) 
norm=stats.shapiro(NTval.dropna())
print(norm)
norm2=stats.shapiro(TBval.dropna())
print(norm2)
norm3=stats.shapiro(TSAval.dropna())
print(norm3)

'''
if Plot :
    fig = plt.figure(figsize= (10, 10))
    ax = fig.add_subplot(111)

    normality_plot, stat = stats.probplot(df['NT'].dropna(), plot= plt, rvalue= True)
    ax.set_title("Probability plot of model residual's", fontsize= 20)
    ax.set

    plt.show() 

    plt.hist(df['NT'].dropna(), bins='auto', histtype='bar', ec='k') 
    plt.xlabel("Residuals")
    plt.ylabel('Frequency')
    plt.show() # Residual =experimental error
'''

# 3.3 homogeneity (check statistically significant difference in their varability)

hmgn=stats.levene(np.squeeze(NTval.dropna()),np.squeeze(TBval.dropna()), np.squeeze(TSAval.dropna())) 
print(hmgn)


#Stage 4 : do the statistic test (2 parameters) depending on parametric test conditions results
##############################

if float(norm[1])<0.05 or float(norm2[1])<0.05 or float(norm3[1])<0.05 or float(hmgn[1])<0.05:
    print ("conditions test failed, let's use a non parametric test (KW)")
    #if conditions NOT checked --> KW test
    kr=stats.kruskal(NTval.dropna(),TBval.dropna(), TSAval.dropna())
    print(kr)
    if float(kr[1])<0.05:
        print("KW test indicates significatives differences between groups values, let's use a post hoc test (Mann Whitney)")
        # Post Hoc test : Mann Whitney test (if KW's pvalue <0.05) 
        mw=stats.mannwhitneyu(NTval.dropna(), TBval.dropna(), alternative='two-sided')
        mw2=stats.mannwhitneyu(NTval.dropna(), TSAval.dropna(), alternative='two-sided')
        mw3=stats.mannwhitneyu(TSAval.dropna(),TBval.dropna(), alternative='two-sided')

        result = {'NT':[float(1.000000),float(mw[1]),float(mw2[1])],'TB':[float(mw[1]),float(1.000000),float(mw3[1])],'TSA':[float(mw2[1]),float(mw3[1]),float(1.000000)]}
        result=pd.DataFrame(result,index=['NT','TB','TSA'])
        print(result)

        if Plot : 
            #Format: diagonal, non-significant, p<0.001, p<0.01, p<0.05
            cmap = ['1', '#c19792',  '#700064','#a00000', '#ec6f6a']
            heatmap_args = {'cmap': cmap, 'linewidths': 0.25, 'linecolor': '0.5', 'clip_on': False, 'square': True, 'cbar_ax_bbox': [0.80, 0.35, 0.04, 0.3]}
            sp.sign_plot(result, **heatmap_args)
            plt.title('Post hoc Mann Whitney test result \n on group comparison of \n mean rise time for images +/=/-',loc='left')
            plt.show()

    else :
        print("KW test indicates NO significative difference between groups values")

else :
    print ("conditions checked, let's use a parametric test (ANOVA one Way)")
    #if conditions checked --> ANOVA test
    pval=Anova(NTval,TBval,TSAval,df2)
    if float(pval) < 0.05 :
        print("Anova test indicates significatives differences between groups values, let's use a post hoc test (Tukey)")
        # Post Hoc test : tukey (if ANOVA's pvalue <0.05) 
        print(stats.tukey_hsd(np.squeeze(NTval.dropna()),np.squeeze(TBval.dropna()), np.squeeze(TSAval.dropna())))
    else : 
        print("Anova test indicates NO significative difference between groups values")



