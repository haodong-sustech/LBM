#include<stdio.h>
#include<stdlib.h>
#include<math.h>
int main()
{
    int nx=201, ny=21;       //x,y方向的格子数
    int i,j,a,a1,ts,ia,ja;
    float **f[9], **ft[9];   //两个指针,分布函数及其时间步更新后的临时存储
    int isn[nx][ny];         //索引
    float wt[9], ux, uy, rho, term1, term2, u2;
    int ex[9], ey[9], kb[9];  //ex, ey为离散速度分量，kb用于处理边界条件反弹后的方向索引
    float source[9];
    int time=20000;           //时间步
    float dpdx=1.0e-5;        //pressure gradient
    float tau=0.8;            //relaxation time
    float rho0=1.0;
    float visc, H;           //粘度系数；v=cs*cs*(tau-0.5)
    float feq[9];
    float rho_av;             //平均密度
    float ux_exact[ny], ux_exact_1[ny];       //用于储存每个y位置的精确（解析解）速度分布，用于验证和比较模拟结果的准确性
    int kpor;                 //用于计算流体格点的总数
    float delta=0.5;
    float y,y2;
    float c1,c2;
    float l2;
    FILE *out, *data;

    visc = (tau-0.5)/3.0;
// H = ny-3.0+3.0*delta;  

    H = ny-1-2.0*delta;

    ex[0] = 0; ey[0] = 0;
    ex[1] = 1; ey[1] = 0;
    ex[2] = 0; ey[2] = 1;
    ex[3] = -1; ey[3] = 0;
    ex[4] = 0; ey[4] = -1;
    ex[5] = 1; ey[5] = 1;
    ex[6] = -1; ey[6] = 1;
    ex[7] = -1; ey[7] = -1;
    ex[8] = 1; ey[8] = -1;

    for (a=0; a<9; a++){               //内存分配
        f[a] = (float**)malloc(nx*sizeof(float*));
            ft[a] = (float**)malloc(nx*sizeof(float*));
        for (i=0; i<nx; i++){
            f[a][i]=(float*)malloc(ny*sizeof(float));
            ft[a][i]=(float*)malloc(ny*sizeof(float));
        }
    }

    for (a=0;a<9;a++){
        if (a==0){wt[a]=4.0/9.0;}      //权重因子
        if (a>=1 & a<=4) {wt[a]=1.0/9.0;}
        if (a>=5 & a<=8) {wt[a]=1.0/36.0;}
    }

        for(a=0;a<9;a++){
            for (a1=a; a1<9; a1++){
                if ((ex[a]+ex[a1])==0 && (ey[a]+ey[a1])==0){
                    kb[a] = a1;
                    kb[a1] = a;   //bounce back 边界反弹
                }
            }
        }
    for (a=0;a<9;a++){
        printf("a = %d, kb = %d\n", a, kb[a]);
    }

// initialize the distribution function
    for(i=0; i<nx; i++){
        for(j=0; j<ny; j++){
            isn[i][j]=0;            //初始化为流体（数值为0），indicator function =0 for fluid and 1 for solid
            if(j==0|j==ny-1){isn[i][j]=1;}   //逻辑运算或是‘|’, 非是“||”
            for(a=0; a<9; a++){
                f[a][i][j]=wt[a]*rho0;     //corresponds to zero velocity
            }
        }
    }

    for (ts=1;ts<=time;ts++) {      //time loop    
        rho_av=0.0; kpor = 0;
        for(i=0;i<nx;i++){
            for(j=1;j<ny-1;j++){    //避免壁面上的格点
                ux=0.0; uy=0.0; rho=0.0;  //初始化
                for(a=0; a<9; a++){
                    rho += f[a][i][j];  //计算密度
                    ux += f[a][i][j]*ex[a];  //计算速度
                    uy += f[a][i][j]*ey[a];
                }
                rho_av += rho;     
                kpor += 1;              //累计参与计算的格子数量
                ux += dpdx/2.0;         //梯度压差引起的速度变化
                ux /= rho; uy /= rho;   //速度归一化，ux,uy是通过分布函数的速度分量累加得到的，因此实际上是动量的形式
                u2 = ux*ux + uy*uy;
                //ux += tau*dpdx;
                for(a=0;a<9;a++){
                    term1 = ux*ex[a] + uy*ey[a];
                    term2 = term1*term1;
                    source[a] = (1.0-0.5/tau)*wt[a]*(3.0*(ex[a]-ux)+9.0*(ex[a]*ux+ey[a]*uy)*ex[a])*dpdx;
                    feq[a] = wt[a]*rho*(1.0+3.0*term1+4.5*term2-1.5*u2);
                    ft[a][i][j] = f[a][i][j]-(f[a][i][j]-feq[a])/tau+source[a];  // collision step
                }
            }
        }
        
        if(fmod(ts,100)==0){
            printf("ts=%d\t rho_av=%12.8f\n", ts, rho_av/kpor);
        }

        //streaming of post-collison particle distributions
        for(i=0;i<nx;i++){
            for(j=1;j<ny-1;j++){
                for(a=0;a<9;a++){
                    ia=i+ex[a];
                    ja=j+ey[a];
                    if(ia<0){ia=nx-1;}   //pushing operator
                    if(ia>nx-1){ia=0;}
                    f[a][ia][ja]=ft[a][i][j];  //a location of (i,j) to a new location (ia,ja)
                }
            }
        }

        //boundary confitions, especially enforcement of bounceback
        for(i=0;i<nx;i++){
            for(j=0;j<ny;j++){
                if(isn[i][j]==0){         //保证计算区域为流体
                    for(a=0;a<9;a++){
                        ia=i-ex[a];       //计算邻近格点的位置
                        ja=j-ey[a];
                        if(ia<0){ia=nx-1;}   //周期性边界条件，左右边界条件
                        if(ia>nx-1){ia=0;}
                        
                        if(isn[ia][ja]==1){
                            f[a][i][j]= f[kb[a]][ia][ja];    //实施反弹性边界条件，上下边界条件
                        }
                    }
                }
            }
        }
    }   

// c1 = 0.5*dpdx*(ny-1)/visc;
// c2 = 0.125*dpdx/visc + 0.25*dpdx*(ny-1)/visc;
    for(j=0;j<ny;j++){         //解析解
        y=(j-delta)/H;
        y2=y*y;
        ux_exact[j]=0.5*dpdx*H*H*(y-y2)/visc;     //等价于ux_exact_1[j]=(H*H/2/visc)*dpdx*((j-delta)/H)*(1-(j-delta)/H);
//      ux_exact[j]=0.5*dpdx*(j-1.0+delta)*(H-(j-(1.0-delta)))/visc;
        if(j==0 | j==ny-1) {ux_exact[j] = 0;}
    
    }

    out = fopen("ux_prof.dat", "w"); //打开文件.dat，以写入的形式
    i=nx-1;
    l2=0.0; c1=0.0; c2=0.0;
    for(j=0;j<ny;j++){
        rho=0.0; ux=0.0; uy=0.0;
        if(isn[i][j]==0){
            for(a=0;a<9;a++){
                rho += f[a][i][j];
                ux += f[a][i][j]*ex[a];
                uy += f[a][i][j]*ey[a];
            }
            ux /= rho; uy /= rho;
        }
        c1 += ux_exact[j]*ux_exact[j];
        c2 += (ux_exact[j]-ux)*(ux_exact[j]-ux);  //误差
        fprintf(out, "%d %d %12.8f %12.8f %12.8f %12.8f\n", nx-1, j, ux, ux_exact[j], uy, rho);  //输出文件的内容，每一列对应的顺序
    }
    printf("c1 = %12.8e c2 = %12.8e\n", c1, c2);
    l2 = pow((c2/c1),0.5);
    printf("l2 = %12.8e\n", 12);
    fclose(out);
}