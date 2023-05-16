/*********************************************
 * OPL 12.8.0.0 Model
 * Author: admin
 * Creation Date: May 13, 2023 at 7:45:22 PM
 *********************************************/
/*********************************************
 * OPL 12.8.0.0 Model
 * Author: admin
 * Creation Date: May 1, 2023 at 7:15:44 PM
 *********************************************/
int r = 5;
int p = 3;
range region = 1..r;
range part = 1..p;
float T[part] = [0.29, 0.29, 0.29];
int n[part] = [10,20,70];
float d[part][region] = [[39.96, 16.49, 8.20, 7.50, 2.47],
		         [0.32, 2.14, 4.16, 5.16, 8.13],
		         [1.12, 2.58, 2.08, 2.56, 3.35]];
float h = 0.2;
float s[part][region] = [[6.98, 6.50, 5.30, 3.50, 4.50],
		         [3.20, 6.20, 6.40, 6.80, 3.55],
		         [2.00, 1.40, 2.40, 3.80, 4.00]];
int l = 5;
int IT = 10;
float CSL = 1.64485362695147;
dvar float+ std[part];
dvar float+ IC[part];
dvar float+ I[part];
dvar float+ ss[part];

minimize sum(i in part, j in region)(365*d[i][j]*n[i]*(I[i]+ T[i]));

subject to{
forall(i in part){
//constraint 1: calculate standard during leadtime
std[i] == sqrt(sum(j in region)((l + IT)*(s[i][j])^2));
}
//constraint 2: calculate inventory cycle 
forall(i in part){
IC[i] == sum(j in region)((d[i][j]*IT)/2);
//constraint 3: calculate safety stock
ss[i] == CSL*std[i];
//constraint 4: calculate inventory cost per unit
I[i] == ((ss[i] + IC[i])/sum(j in region)d[i][j])*h;
}
}
 