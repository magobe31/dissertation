
%%%% Resolution of Subproblem 2 - Optimization of a single reflection
%%%% coefficient with all the other variables being fixed 

function [Am,Bm,Final_matrix,Am1,Bm1,NZ_eigenvalues,optimal_solution_1,optimal_capacity_1] =  compute_Am_Bm(M,R,T,H,optR,optQ,Nr,Nt,sigma_lineal)

[U1,D1,V1] = eig(optQ);

H_prima = H*U1*D1^1/2;

B = zeros(Nr,Nt);
C = zeros(Nr,Nt);

 

Am1 = [];
Bm1 = [];
Final_matrix = [];
S111 = [];  
U111 = [];
Vm1 = [];
Vm_inv_1 = [];
vm_vec = [];
vmT_fin = [];
optimal_solution_1 = [];
optimal_capacity_1 = [];

  

    for l = 1:M   
    
         for k = 1:M
               if k ~= l
                    B1 = optR(k)*R(:,k)*((D1^1/2*V1*T(k,:)')');
                    B = B + B1;
                    Am = eye(Nr) + 1/sigma_lineal*(H_prima + B)*(H_prima + B)' + 1/sigma_lineal * R(:,l) * ((D1^1/2*V1*T(k,:)')') * (D1^1/2*V1*T(k,:)') * R(:,l)';
      
                    C1 = (D1^1/2*V1*T(k,:)')*R(:,k)'*conj(optR(k));
                    C = C + C1;
                    Bm = 1/sigma_lineal * R(:,l) * ((D1^1/2*V1*T(l,:)')') * (H_prima' + C);
               end
         end    
        
        Am1 = [Am1 Am];
        Bm1 = [Bm1 Bm];
        
        Final_matrix_1 = (Am1(:, 1+(l-1)*Nr:Nr*l))\Bm1(:, 1+(l-1)*Nr:Nr*l); 
        Final_matrix = [Final_matrix Final_matrix_1];
    
        [U11,S11,V11] = eig(Final_matrix(:, 1+(l-1)*Nr:Nr*l));
        S111 = [S111  S11(1)];
        U111 = [U111 U11];
        NZ_eigenvalues = S111;
        Vm = U111(:, 1+(l-1)*Nr:Nr*l)'*Am1(:, 1+(l-1)*Nr:Nr*l)*U111(:, 1+(l-1)*Nr:Nr*l);
        Vm1 = [Vm1 Vm];
        
        Vm_inv = inv(Vm1(:, 1+(l-1)*Nr:Nr*l));
        Vm_inv_1 = [Vm_inv_1 Vm_inv];
        vm = Vm_inv_1(:,1+(l-1)*Nr);
        vm_vec = [vm_vec vm];
        vm1 = vm_vec(1,:);   %%% first elements of vm_vec
        vmT = Vm1(1,:);
        vmT1 = vmT(1+(l-1)*Nr);
        vmT_fin = [vmT_fin vmT1];  %%% first element of vmT1
       
        
        if trace(Final_matrix(:, 1+(l-1)*Nr:Nr*l)) ~= 0 
            
            optimal_solution = exp(-1i*angle(NZ_eigenvalues(l)));
            optimal_capacity = log2(1+abs(NZ_eigenvalues(l))^2*(1-vmT_fin(l)*vm1(l))+2*abs(NZ_eigenvalues(l)))+log2(det(Am1(:, 1+(l-1)*Nr:Nr*l)));
            
        else
            optimal_solution = 1;  
            optimal_capacity = log2(det(Am1(:, 1+(l-1)*Nr:Nr*l)-Bm1(:, 1+(l-1)*Nr:Nr*l)'*inv(Am1(:, 1+(l-1)*Nr:Nr*l))*Bm1(:, 1+(l-1)*Nr:Nr*l)));
        end
       
       optimal_solution_1 = [optimal_solution_1 optimal_solution];    
       optimal_capacity_1 = [optimal_capacity_1 optimal_capacity];
       
   end

end