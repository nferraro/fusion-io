function jac = jacobian3(nr, nphi, ntheta, R,Z,nfp, theta,phi,psi)
nargin
  if(nr==1)       
    disp("Error: cannot calculate jacobian with one surface");
  end
  jac_sign = 0;

  % Check if last three arguments are provided
  if nargin < 9
    disp("Last three arguments (theta, phi, psi) are missing.Calculating Jacobian as in fusion-io");
  end

  

if (nargin == 6)  
dphidj =  2.* pi / nphi/nfp;
for i=1:nr
  for j=1:nphi
     for k=1:ntheta
         jac_sign = 0;
        if(i==1)
            dRdi = (R(i+1,j,k) - R(i,j,k));
	        dZdi = (Z(i+1,j,k) - Z(i,j,k));
        elseif (i==nr)
            dRdi = (R(i,j,k) - R(i-1,j,k));
	        dZdi = (Z(i,j,k) - Z(i-1,j,k));
        else
            dRdi = (R(i+1,j,k) - R(i-1,j,k))/2;
	        dZdi = (Z(i+1,j,k) - Z(i-1,j,k))/2;
        end

        if(j==1)
            dRdj = (R(i,j+1,k) - R(i,j,k));
	        dZdj = (Z(i,j+1,k) - Z(i,j,k));
        elseif (j==nphi)
            dRdj = (R(i,j,k) - R(i,j-1,k));
	        dZdj = (Z(i,j,k) - Z(i,j-1,k));
        else
            dRdj = (R(i,j+1,k) - R(i,j-1,k))/2;
	        dZdj = (Z(i,j+1,k) - Z(i,j-1,k))/2;
        end


        if(k==1)
            dRdk = (R(i,j,k+1) - R(i,j,k));
	        dZdk = (Z(i,j,k+1) - Z(i,j,k));
        elseif (k==ntheta)
            dRdk = (R(i,j,k) - R(i,j,k-1));
	        dZdk = (Z(i,j,k) - Z(i,j,k-1));
        else
            dRdk = (R(i,j,k+1) - R(i,j,k-1))/2;
	        dZdk = (Z(i,j,k+1) - Z(i,j,k-1))/2;
        end

	
	jac(i,j,k)=(dRdi*dZdk - dRdk*dZdi)*R(i,j,k)*dphidj;
	
	% // di . (dk x dj)
	% // where i is the radial index
	% //       k is the poloidal index
	% //       j is the toroidal index


	if(jac_sign==0)
        if(jac(i,j,k)>0)
            jac_sign=1;
        else
            jac_sign=-1;
        end
	else 
	  if(jac(i,j,k)*jac_sign <= 0.) 
	   disp("Error: Jacobian flips sign")
       i
       j
       k
      end
    end
	
     end
  end
end
elseif(nargin == 9)
for i=1:nr
 for j=1:nphi
     for k=1:ntheta
jac_sign = 0;
        if(i==1)
            dRdpsi = (R(i+1,j,k) - R(i,j,k))/(psi(i+1,1)-psi(i,1));
	        dZdpsi = (Z(i+1,j,k) - Z(i,j,k))/(psi(i+1,1)-psi(i,1));
        elseif (i==nr)
            dRdpsi = (R(i,j,k) - R(i-1,j,k))/(psi(i,1)-psi(i-1,1));
	        dZdpsi = (Z(i,j,k) - Z(i-1,j,k))/(psi(i,1)-psi(i-1,1));
        else
            dRdpsi = (R(i+1,j,k) - R(i-1,j,k))/(psi(i+1,1)-psi(i-1,1));
	        dZdpsi = (Z(i+1,j,k) - Z(i-1,j,k))/(psi(i+1,1)-psi(i-1,1));
        end

        if(j==1)
            dRdphi = (R(i,j+1,k) - R(i,1,k))/(phi(1,j+1)-phi(1,j));
	        dZdphi = (Z(i,j+1,k) - Z(i,1,k))/(phi(1,j+1)-phi(1,j));
        elseif (j==nphi)
            dRdphi = (R(i,j,k) - R(i,j-1,k))/(phi(1,j)-phi(1,j-1));
	        dZdphi = (Z(i,j,k) - Z(i,j-1,k))/(phi(1,j)-phi(1,j-1));
        else
            dRdphi = (R(i,j+1,k) - R(i,j-1,k))/(phi(1,j+1)-phi(1,j-1));
	        dZdphi = (Z(i,j+1,k) - Z(i,j-1,k))/(phi(1,j+1)-phi(1,j-1));
        end


        if(k==1)
            dRdtheta = (R(i,j,k+1) - R(i,j,k))/(theta(1,k+1)-theta(1,k));
	        dZdtheta = (Z(i,j,k+1) - Z(i,j,k))/(theta(1,k+1)-theta(1,k));
        elseif (k==ntheta)
            dRdtheta = (R(i,j,k) - R(i,j,k-1))/(theta(1,k)-theta(1,k-1));
	        dZdtheta = (Z(i,j,k) - Z(i,j,k-1))/(theta(1,k)-theta(1,k-1));
        else
            dRdtheta = (R(i,j,k+1) - R(i,j,k-1))/(theta(1,k+1)-theta(1,k-1));
	        dZdtheta = (Z(i,j,k+1) - Z(i,j,k-1))/(theta(1,k+1)-theta(1,k-1));
        end


        % jac(i,j,k)=(dRdi*dZdk - dRdk*dZdi)*R(i,j,k)*dphidj;
        % [del psi.(del theta X del phi)]^-1
        % Jac_vmec(i,j,k)=dpsidi*dthetadk - dpsidk*dthetadi)*dphidj;

        jac(i,j,k)=(dRdtheta*dZdpsi-dRdpsi*dZdtheta)*R(i,j,k);
        
	    if(jac_sign==0)
            if(jac(i,j,k)>0)
                jac_sign=1;
            else
                jac_sign=-1;
            end
	    else 
	      if(jac(i,j,k)*jac_sign <= 0.) 
	       disp("Error: Jacobian flips sign")
           i
           j
           k
          end
        end
	
     end
 end
end
end

