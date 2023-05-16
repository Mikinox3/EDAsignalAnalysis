import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
import numpy as np
from bioinfokit.analys import stat
from statsmodels.graphics.factorplots import interaction_plot

#Stage 1 : get the file
##############################
 
aFile = 'C:/xxxx/MPeakAmp2.csv'
df = pd.read_csv(aFile,sep=';', header=0)
print(df)

#Stage 2 : select what you want to compare between groups 
##############################
select = ['Ipos','Ineu','Ineg'] # group comparison for img all valences
#select = ['Fpos','Fneu','Fneg'] # group comparison for negative img 
df2 = df[df['Valence'].isin(select)]
print(df2)

#prepare the data for stats
NTval=df2.loc[df2['Groupe']=='NT',['Value']]
TBval=df2.loc[df2['Groupe']=='TB',['Value']]
TSAval=df2.loc[df2['Groupe']=='TSA',['Value']]

#plot the data
sns.boxplot(x="Valence", y="Value", hue="Groupe", data=df2, palette="Set3") 
plt.title('val : group comparison of images by valence')
plt.show()


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

#Stage 4 : do the statistic test (3 parameters) depending on parametric test conditions results
##############################


#if float(norm[1])<0.05 or float(norm2[1])<0.05 or float(norm3[1])<0.05 or float(hmgn[1])<0.05:
    #print ("conditions test failed, let's use a non parametric test (KW)")
#else :

#print("conditions checked, let's use a parametric test (ANOVA 2 way)")
# ANOVA 2 way test
model = ols('Value ~ C(Valence) + C(Groupe) + C(Valence):C(Groupe)', data=df2).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print(anova_table)

fig = interaction_plot(x=df2['Valence'].to_numpy(), trace=df2['Groupe'].to_numpy(), response=df2['Value'].to_numpy(), colors=['#4c061d','#d17a22', '#b4c292'],markers=['D','^','*'],xlabel='Valence',ylabel='val')
plt.title("Interaction plot :  groupe + valence comparison of the val")
plt.show()

if anova_table["PR(>F)"][0]<0.05 or anova_table["PR(>F)"][1]<0.05 or anova_table["PR(>F)"][2]<0.05:

    print("Anova test indicates significatives differences between groups and valences values, let's use a post hoc test (Tukey)")

    # Post hoc test (Tukey) if p-value<0.05
    res=stat()
    # check valence effect on value
    res.tukey_hsd(df=df2, res_var='Value', xfac_var='Valence', anova_model='Value ~ C(Valence) + C(Groupe) + C(Valence):C(Groupe)')
    print(res.tukey_summary)

    #check group effect on value
    res.tukey_hsd(df=df2, res_var='Value', xfac_var='Groupe', anova_model='Value ~ C(Valence) + C(Groupe) + C(Valence):C(Groupe)')
    print(res.tukey_summary)

    # check interaction effect between groupe and valence
    res.tukey_hsd(df=df2, res_var='Value', xfac_var=['Groupe','Valence'], anova_model='Value ~ C(Valence) + C(Groupe) + C(Valence):C(Groupe)')
    print(res.tukey_summary)
else :
    print("Anova test indicates NO significative difference between groups and valence values")




