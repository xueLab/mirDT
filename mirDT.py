# mirDT is a decision-tree based meta-predictor for miRNA prediction by integrating
#       the predictive results of MiPred, ProMiR and miRPara.
# USAGE: ./mirDT.py
#       After typing in the above-mentioned command, the user will be asked to type 
#             in the name of inpout file.
#       Input file is a tab-delimited file with containing miRNA_ID and its 
#             predictive scores from miPred, ProMiR, and miRPara, respectively.
#

import math

def read_input_file(filename):
    x = open(filename,"rt")
    a = x.readlines()
    input_lst = []
    for i in a:
        if i.split("\t")[2]=="NA":
            input_lst.append(i.split("\t")[0]+"\t"+i.split("\t")[1]+"\t"+"-15"+"\t"+i.split("\t")[3])
        elif i.split("\t")[2]=="0.00000000000":
            input_lst.append(i.split("\t")[0]+"\t"+i.split("\t")[1]+"\t"+"-12"+"\t"+i.split("\t")[3])
        else:
            #print float(i.split("\t")[2])
            a=math.log(float(i.split("\t")[2]))
            #print float(i.split("\t")[2])
            input_lst.append(i.split("\t")[0]+"\t"+i.split("\t")[1]+"\t"+str(a)+"\t"+i.split("\t")[3])
    return input_lst

def decision_tree_first(input_lst,Tp1,Tp2,Tp3,Tn1,Tn2,Tn3):
    pos_lst = []
    neg_lst = []
    non_lst = []
    Tp = [Tp1,Tp2,Tp3]
    Tn = [Tn1,Tn2,Tn3]
    for i in input_lst:
        if float(i.split("\t")[1].strip())>=Tp1 and float(i.split("\t")[2].strip())>=Tp2 and float(i.split("\t")[3].strip())>=Tp3:
            pos_lst.append(i)
        elif float(i.split("\t")[1].strip())>=Tp1 and float(i.split("\t")[2].strip())>=Tp2:
            pos_lst.append(i)
        elif float(i.split("\t")[1].strip())>=Tp1 and float(i.split("\t")[3].strip())>=Tp2:
            pos_lst.append(i)
        elif float(i.split("\t")[2].strip())>=Tp1 and float(i.split("\t")[3].strip())>=Tp2:
            pos_lst.append(i)
        elif float(i.split("\t")[1].strip())<=Tn1 and float(i.split("\t")[2].strip())<=Tn2 and float(i.split("\t")[3].strip())<=Tn3:
            neg_lst.append(i)
        elif float(i.split("\t")[1].strip())<=Tn1 and float(i.split("\t")[2].strip())<=Tn2:
            neg_lst.append(i)
        elif float(i.split("\t")[1].strip())<=Tn1 and float(i.split("\t")[3].strip())<=Tn2:
            neg_lst.append(i)
        elif float(i.split("\t")[2].strip())<=Tn1 and float(i.split("\t")[3].strip())<=Tn2:
            neg_lst.append(i)
        elif float(i.split("\t")[1].strip())>=Tp1 and float(i.split("\t")[2].strip())<=Tn2:
            dist_pos = (float(i.split("\t")[1].strip())-Tp[0])**2+(float(i.split("\t")[2].strip())-Tp[1])**2+(float(i.split("\t")[3].strip())-Tp[2])**2
            dist_neg = (float(i.split("\t")[1].strip())-Tn[0])**2+(float(i.split("\t")[2].strip())-Tn[1])**2+(float(i.split("\t")[3].strip())-Tn[2])**2
            if dist_pos>dist_neg:
                pos_lst.append(i)
            elif dist_pos<dist_neg:
                neg_lst.append(i)
        elif float(i.split("\t")[1].strip())>=Tp1 and float(i.split("\t")[3].strip())<=Tn3:
            dist_pos = (float(i.split("\t")[1].strip())-Tp[0])**2+(float(i.split("\t")[2].strip())-Tp[1])**2+(float(i.split("\t")[3].strip())-Tp[2])**2
            dist_neg = (float(i.split("\t")[1].strip())-Tn[0])**2+(float(i.split("\t")[2].strip())-Tn[1])**2+(float(i.split("\t")[3].strip())-Tn[2])**2
            if dist_pos>dist_neg:
                pos_lst.append(i)
            elif dist_pos<dist_neg:
                neg_lst.append(i)
        elif float(i.split("\t")[2].strip())>=Tp2 and float(i.split("\t")[1].strip())<=Tn1:
            dist_pos = (float(i.split("\t")[1].strip())-Tp[0])**2+(float(i.split("\t")[2].strip())-Tp[1])**2+(float(i.split("\t")[3].strip())-Tp[2])**2
            dist_neg = (float(i.split("\t")[1].strip())-Tn[0])**2+(float(i.split("\t")[2].strip())-Tn[1])**2+(float(i.split("\t")[3].strip())-Tn[2])**2
            if dist_pos>dist_neg:
                pos_lst.append(i)
            elif dist_pos<dist_neg:
                neg_lst.append(i)
        elif float(i.split("\t")[2].strip())>=Tp2 and float(i.split("\t")[3].strip())<=Tn3:
            dist_pos = (float(i.split("\t")[1].strip())-Tp[0])**2+(float(i.split("\t")[2].strip())-Tp[1])**2+(float(i.split("\t")[3].strip())-Tp[2])**2
            dist_neg = (float(i.split("\t")[1].strip())-Tn[0])**2+(float(i.split("\t")[2].strip())-Tn[1])**2+(float(i.split("\t")[3].strip())-Tn[2])**2
            if dist_pos>dist_neg:
                pos_lst.append(i)
            elif dist_pos<dist_neg:
                neg_lst.append(i)
        elif float(i.split("\t")[3].strip())>=Tp3 and float(i.split("\t")[1].strip())<=Tn1:
            dist_pos = (float(i.split("\t")[1].strip())-Tp[0])**2+(float(i.split("\t")[2].strip())-Tp[1])**2+(float(i.split("\t")[3].strip())-Tp[2])**2
            dist_neg = (float(i.split("\t")[1].strip())-Tn[0])**2+(float(i.split("\t")[2].strip())-Tn[1])**2+(float(i.split("\t")[3].strip())-Tn[2])**2
            if dist_pos>dist_neg:
                pos_lst.append(i)
            elif dist_pos<dist_neg:
                neg_lst.append(i)
        elif float(i.split("\t")[3].strip())>=Tp3 and float(i.split("\t")[2].strip())<=Tn2:
            dist_pos = (float(i.split("\t")[1].strip())-Tp[0])**2+(float(i.split("\t")[2].strip())-Tp[1])**2+(float(i.split("\t")[3].strip())-Tp[2])**2
            dist_neg = (float(i.split("\t")[1].strip())-Tn[0])**2+(float(i.split("\t")[2].strip())-Tn[1])**2+(float(i.split("\t")[3].strip())-Tn[2])**2
            if dist_pos>dist_neg:
                pos_lst.append(i)
            elif dist_pos<dist_neg:
                neg_lst.append(i)
        elif float(i.split("\t")[1].strip())>=float(Tp1):
            pos_lst.append(i)
        elif float(i.split("\t")[2].strip())>=float(Tp2):
            pos_lst.append(i)
        elif float(i.split("\t")[3].strip())>=float(Tp3):
            pos_lst.append(i)
        elif float(i.split("\t")[1].strip())<=float(Tn1):
            neg_lst.append(i)
        elif float(i.split("\t")[2].strip())<=float(Tn2):
            neg_lst.append(i)
        elif float(i.split("\t")[3].strip())<=float(Tn3):
            neg_lst.append(i)
        else:
            non_lst.append(i)
    return pos_lst,neg_lst,non_lst

def decide_second_neg(non_lst,Tn1,Tn2,Tn3):
    pos_lst_2 = []
    neg_lst = []
    for i in non_lst:
        if float(i.split("\t")[1].strip())<=float(Tn1) and float(i.split("\t")[3].strip())<=float(Tn3):
            neg_lst.append(i)
        elif float(i.split("\t")[1].strip())<=float(Tn1) and float(i.split("\t")[2].strip())<=float(Tn2):
            neg_lst.append(i)
        elif float(i.split("\t")[3].strip())<=float(Tn3) and float(i.split("\t")[2].strip())<=float(Tn2):
            neg_lst.append(i)
        else:
            pos_lst_2.append(i)
    neg_lst_2 = []
    [neg_lst_2.append(j) for j in neg_lst if j not in neg_lst_2]
    return pos_lst_2, neg_lst_2

def decide_second_pos(non_lst,Tp1,Tp2,Tp3):
    pos_lst = []
    neg_lst_2 = []
    for i in non_lst:
        if float(i.split("\t")[1].strip())>=float(Tp1) and float(i.split("\t")[3].strip())>=float(Tp3):
            pos_lst.append(i)
        elif float(i.split("\t")[1].strip())>=float(Tp1) and float(i.split("\t")[2].strip())>=float(Tp2):
            pos_lst.append(i)
        elif float(i.split("\t")[3].strip())>=float(Tp3) and float(i.split("\t")[2].strip())>=float(Tp2):
            pos_lst.append(i)
        else:
            neg_lst_2.append(i)
    pos_lst_2 = []
    [pos_lst_2.append(j) for j in pos_lst if j not in pos_lst_2]
    return pos_lst_2, neg_lst_2

def main(filename):
    input_lst = read_input_file(filename)
    TP_lst = open("TP_lst.txt","wt")
    TN_lst = open("TN_lst.txt","wt")
    pos_lst,neg_lst,non_lst = decision_tree_first(input_lst,58.40,-2.36,0.967,52.12,-8.06,0.796)
    if len(pos_lst)!=0:
        TP = pos_lst
    elif len(pos_lst)==0:
        TP_lst_2, nonFN_lst_2 = decide_second_pos(non_lst,52.12,-2.82,0.923)
        TP = pos_lst + TP_lst_2
    if len(neg_lst)!=0:
        TN = neg_lst
    elif len(neg_lst)==0:
        TN_lst_2, nonFP_lst_2 = decide_second_neg(non_lst,60.4,-6.25,0.875)
        TN = neg_lst + TN_lst_2
    for i in TP:
        TP_lst.write(i)
    for j in TN:
        TN_lst.write(j)
    return TP, TN

filename = raw_input("Enter a Filename:")
main(filename)


