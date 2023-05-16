	/*********************************************
	 * OPL 12.8.0.0 Model
	 * Author: admin
	 * Creation Date: May 11, 2023 at 12:39:14 AM
	 *********************************************/
	int r = 5;
	int p = 3;
	range region = 1..r;
	range part = 1..p;
	float T1[part] = [0.29, 0.29, 0.29];
	float T2[region] = [0.19, 0.19, 0.19, 0.19, 0.19];
	int n[part] = [10, 20, 70];
	float d[part][region] = [
	    [39.96, 16.49, 8.20, 7.50, 2.47],
	    [0.32, 2.14, 4.16, 5.16, 8.13],
	    [1.12, 2.58, 2.08, 2.56, 3.35]
	];
	float h = 0.2;
	float s[part][region] = [
	    [6.98, 6.50, 5.30, 3.50, 4.50],
	    [3.20, 6.20, 6.40, 6.80, 3.55],
	    [2.00, 1.40, 2.40, 3.80, 4.00]
	];
	int numdc = 5;
	int lc = 5;
	int ld = 7;
	int IT = 10;
	float CSL = 1.64485362695147;
	// Centralization
	dvar float+ x1[part];
	dvar boolean y1[part];
	dvar float+ std1[part];
	dvar float+ IC1[part];
	dvar float+ I1[part];
	dvar float+ ss1[part];
	
	// Decentralization
	dvar boolean y2[part][region];
	dvar float+ x2[part][region];
	dvar float+ std2[part][region];
	dvar float+ IC2[part][region];
	dvar float+ I2[part][region];
	dvar float+ ss2[part][region];
	
	//decide to centralization or decentralization 
	dvar boolean centralized[part];
	dvar boolean decentralized[part];
	
	// Expressions for inventory management
	dexpr float cd = sum(i in part, j in region)(decentralized[i]*d[i][j]*n[i]*365*(I2[i][j]+T2[j]));
	dexpr float cc = sum(i in part, j in region)(centralized[i]*d[i][j]*n[i]*365*(I1[i]+T1[i]));
	
	// Objective function
	minimize cc + cd;
	
	subject to {
	    // Constraint 1: At least one part must be centralized
	    sum(i in part) centralized[i] >= 1;
	    // Constraint 2: A part can either be centralized or decentralized, but not both
	    forall(i in part) {
    		centralized[i] + decentralized[i] == 1; 
		}
		forall(j in region){
		  sum(i in part) y2[i][j] >= 1;
  		}		  
	    // Constraint 3: If part i is centralized, set decentralization variables to 0 
		forall(i in part) {
    		decentralized[i] <= 1 - centralized[i];
    	}    		
    	forall(i in part,j in region) {
        	y2[i][j] <= 1 - centralized[i];
    	}
		// Constraint 4: If part i is decentralized, set centralization variables to 0
		forall(i in part) {
    		centralized[i] <= 1 - decentralized[i]; 
    		y1[i] <= 1 - decentralized[i];
		} 
	    // Constraint 5: The total number of data centers must be at least numdc
	   sum(i in part) y1[i] + sum(i in part, j in region) y2[i][j] >= numdc;
	
	    // Constraint 6: Compute the centralized demand for each part
	    forall(i in part, j in region) {
        	x2[i][j] == (1 - centralized[i]) * d[i][j];
	    }
	
	    // Constraint 7: Compute the decentralized demand for each part and region
	    forall(i in part, j in region) {
        	x1[i] == (1 - decentralized[i]) * sum(j in region) d[i][j];
	    }
	    // DECENTRALIZATION
	    // Constraint 8: Calculate standard deviation during lead time
	    forall(i in part, j in region) {
	        std2[i][j] == sqrt((ld + IT) * (s[i][j]) ^ 2);
	    }
	
	    // Constraint 9: Calculate inventory cycle
	    forall(i in part, j in region) {
	        IC2[i][j] == (x2[i][j] * IT) / 2;
	        
	        // Constraint 9: Calculate safety stock
	        ss2[i][j] == CSL * std2[i][j];
	        
	        // Constraint 10: Calculate inventory cost per unit
	        I2[i][j] == ((ss2[i][j] + IC2[i][j]) / d[i][j]) * h;
	    }
	
	    // CENTRALIZATION
	    // Constraint 11: Calculate standard deviation during lead time
	    forall(i in part) {
	        std1[i] == sqrt(sum(j in region) ((lc + IT) * (s[i][j]) ^ 2));
	    }
	
	    // Constraint 12: Calculate inventory cycle
	    forall(i in part) {
	        IC1[i] == (x1[i] * IT) / 2;
	
	        // Constraint 13: Calculate safety stock
	        ss1[i] == CSL * std1[i];
	
	        // Constraint 14: Calculate inventory cost per unit
	        I1[i] == ((ss1[i] + IC1[i]) / sum(j in region) d[i][j]) * h;
	    } 
}
