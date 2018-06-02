namespace KM {
		int lenx,leny;
		int w[MAX_NODE][MAX_NODE]; 
		int slack[MAX_NODE],lx[MAX_NODE],ly[MAX_NODE],maty[MAX_NODE];  
		bool vx[MAX_NODE],vy[MAX_NODE]; //S集合、Y集合   
		 
		bool search(int u) {  
		    int i,t;  
		    vx[u]=1;  
		    for(i=0;i<leny;++i)  
		        if(!vy[i]) {  
		            t=lx[u]+ly[i]-w[u][i];  
		            if (t==0) {  
		                vy[i]=1;  
		                if(maty[i]==-1||search(maty[i])){  
		                    maty[i]=u;  
		                    return 1;  
		                }  
		            }  
		            else if(slack[i]>t)  
		                slack[i]=t;  
		        }  
		    return 0;  
		}  
		int get() {
		    int i,j,ans=0;  
		    for(i=0;i<lenx;++i)  
		        for(lx[i]=-INF,j=0;j<leny;++j)  
		            lx[i]=max(lx[i],w[i][j]);  
		    memset(maty,-1,sizeof(maty));  
		    memset(ly,0,sizeof(ly));  
		    for(i=0;i<lenx;++i) { //找增广路 
		        for(j=0;j<leny;++j)  
		            slack[j]=INF;  
		        while(1) {  
		            memset(vx,0,sizeof(vx));  
		            memset(vy,0,sizeof(vy));  
		            if(search(i))//找到i对应的增广路，不再找  
		                break;  
		            //没找到增广路，修正  
		            int d=INF;  
		            for(j=0;j<leny;++j)  
		                if(!vy[j]&&d>slack[j])  
		                    d=slack[j];  
		            for(j=0;j<lenx;++j)  
		                if(vx[j])  
		                    lx[j]-=d;  
		            for(j=0;j<leny;++j)  
		                if(vy[j])  
		                    ly[j]+=d;  
		        }  
		    }  
		    for(i=0;i<leny;++i)  
		        if(maty[i]!=-1)  
		            ans+=w[maty[i]][i];  
		    return ans;  
		}

		void init(const answer& ans) {
			lenx = 0;
			leny = 0;
			vector<int> target_not_match_list;

			for (int i = 0; i < ans.target_map.size(); ++i) {
				if (ans.target_map[i] == 0) {
					leny ++;
					target_not_match_list.push_back(i);
				}
			}

			for (int i = 0; i < ans.match.size(); ++i) {
				 if (ans.match[i] == 0) {
				 		lenx ++;
				 		for (int j = 0; j < target_not_match_list.size(); ++j) {
				 			w[lenx][j + 1] = ans.calc_edit_cost(i, target_not_match_list[j]);
				 		}
				 }
			}
		}
};