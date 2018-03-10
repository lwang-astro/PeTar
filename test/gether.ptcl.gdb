#egrep '[xyz] =' check|sed -e 's/.*[xyz] =\(.*\)/\1/' -e 's/,//g'|awk 'BEGIN{n=0;i=0;printf("x%d=np.array([",i);i++;} {if(n<5) {printf("%s,",$1);n++;} else {n=0;printf("%s])\nx%d=np.array([", $1,i);i++};}'

egrep 'id =' check|sed -e 's/.*id =\(.*\)/\1/' -e 's/,//g'|awk '{print $1}' >id
egrep 'x =' check|sed -e 's/.*x =\(.*\)/\1/' -e 's/,//g'|awk '{if(NR%2!=0) print $1}' >x
egrep 'y =' check|sed -e 's/.*y =\(.*\)/\1/' -e 's/,//g'|awk '{if(NR%2!=0) print $1}' >y
egrep 'z =' check|sed -e 's/.*z =\(.*\)/\1/' -e 's/,//g'|awk '{if(NR%2!=0) print $1}' >z
egrep 'r_search =' check|sed -e 's/.*r_search =\(.*\)/\1/' -e 's/,//g'|awk '{print $1}' >rs
paste id x y z rs |awk '{printf("x[%d]=np.array([%s,%s,%s,%s,%s])\n",NR-1,$1,$2,$3,$4,$5)}' > r
