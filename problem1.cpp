#include<iostream>
#include<fstream>
#include<vector> 
#include<cmath>

int main() {
    //grid size = resolution
    int Nx = 1024; 
    int Ny = 1024;
    double dx = 1.0/(Nx); //dx and dy are the same 
    std::vector<double>psi(Nx*Nx , 0.0);
    std::vector<std::vector<double>> pixels(Ny, std::vector<double> (Nx));
    double omega = 2.0- 1.0/(Nx); //starting omega value 
    double tol = 1e-9; //tolerance

    //boundary conditions
    for(int j = 1; j < Ny-1; j++ ){
        int n = j*Nx;
        psi[n] = psi[n+1]; 
        int i = Nx-1; 
        n = i + j*Nx; 
        psi[n] = psi[n-1];
    }
    for(int i = 1; i < Nx-1; i++){
        int n = i;
        psi[n] = psi[n+Nx];
        int j = Ny-1;
        n = i + j*Nx;
        psi[n] = psi[n-Nx];
    }
    //defining rho 
    std::vector<double>rho(Nx*Nx); 
    for(int j =1; j <Ny-1; j++){
        for(int i = 1; i < Nx-1; i++){
            int n = i+j*Nx;
            double x = (i+1)*dx;
            double y = (j+1)*dx;
            if(0.3<x && x < 0.7 && 0.65 < y && y< 0.7){
                rho[n] = 1.0; 
            } 
            else if(0.3 < x && x < 0.7 && 0.3 < y && y< 0.35){
                rho[n] = -1.0;
            }
            else{
                rho[n] = 0.0;
        }
    }
}
    
    std::vector<double>resid(Nx*Nx, 0.0);
    double sum = 0; 
    for(int k = 0; k < 10000; k++){
        sum = 0; 
        for(int i = 1; i < Nx-1; i++){
            for(int j = 1; j < Ny-1; j++){
                int n = i+j*Nx; 
                psi[n] = psi[n] + 0.25*omega*(psi[n+1] + psi[n-1]+ psi[n+Nx] + psi[n-Nx] + rho[n]*dx*dx - 4* psi[n]); 
                sum = sum + (psi[n+1] + psi[n-1]+ psi[n+Nx] + psi[n-Nx] + rho[n]*dx*dx - 4* psi[n])*(psi[n+1] + psi[n-1]+ psi[n+Nx] + psi[n-Nx] + rho[n]*dx*dx - 4* psi[n]); 
                pixels[j][i] = psi[n];
            }
        }
        //enforce boundary conditions
        for(int j = 1; j < Ny-1; j++ ){
            int n = j*Nx;
            psi[n] = psi[n+1]; 
            int i = Nx-1; 
            n = i + j*Nx; 
            psi[n] = psi[n-1];
        }
        for(int i = 1; i < Nx-1; i++){
            int n = i;
            psi[n] = psi[n+Nx];
            int j = Ny-1;
            n = i + j*Nx;
            psi[n] = psi[n-Nx];
        }
        resid[k] = std::sqrt(sum);
        std::cout << k << " " << resid[k] << std::endl;
        if(resid[k]< tol){
            break;
        }
    }
    std::ofstream out("psi.csv");
    for(int j =0; j<Ny;j++){
          out<<pixels[j][0];
          for(int i = 1; i<Nx;i++){
               out<<","<<pixels[j][i];
          }
          out<<std::endl;
     }
}