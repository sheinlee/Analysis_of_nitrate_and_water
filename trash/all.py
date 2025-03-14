LIST_OF_ELEMENT = ['La','Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu']
for element in LIST_OF_ELEMENT:
    with open('all_BLPHEN_OK.txt','a') as f:
        f.close()
    fi = open(element+'_BLPHEN_OK.txt','r')
    txt = fi.readlines()
    for w in txt:
        with open('all_BLPHEN_OK.txt','a') as f:
            f.write(w)
            f.close()
    fi.close()
