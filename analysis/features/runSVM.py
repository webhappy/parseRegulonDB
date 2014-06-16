import csv,pandas as pp
from sklearn import svm
from sklearn.cross_validation import KFold

df=pp.read_csv('features.csv',index_col=False)
xVals=df[[col for col in df.columns if col != 'Target' and col != 'GeneName' and col !='entireSeq']]
yVals=df['Target']

clf=svm.SVC()


outFile=csv.writer(open('svmResults.csv','wb'))
kf = KFold(len(yVals),n_folds=len(yVals),indices=False)
outFile.writerow(['Gene','Actual','Predict'])
for train,test in kf:
    xTrain=xVals[train]
    xTest=xVals[test]
    yTrain=yVals[train]
    yTest=df[test]
    clf=clf.fit(xTrain,yTrain)
    fitted=clf.predict(xTest)[0]
    actual=yTest['Target'].values[0]
    #print "Predict %f when actually is %f, gene %s at dist %i"%(fitted,actual,yTest['GeneName'].values[0],yTest['DistanceFromATG'].values[0])
    outFile.writerow([yTest['GeneName'].values[0],actual,fitted])
