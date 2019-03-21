%Algorithm is as defined in the book 'Numerical Optimization' p.60
function [alpha_star] = zoom(x, grad, A, b, c1, c2, alpha_lo, alpha_hi)
  
	current_loss = 0.5*x'*A*x - b'*x;
	p = -1.0*grad;
	
	phi_prime_zero = grad'*p;
	
  	epsilon = 0.0001;

	alpha_star = 0.5 * (alpha_hi+alpha_lo);		%initialize alpha_star
	
	while ( abs(alpha_hi-alpha_lo)>epsilon )  
    
		alpha = 0.5 * (alpha_hi+alpha_lo); 
		cost_alpha = 0.5*(x + alpha*p)'*A*(x + alpha*p) - b'*(x + alpha*p);
		cost_alphalo = 0.5*(x + alpha_lo*p)'*A*(x + alpha_lo*p) - b'*(x + alpha_lo*p);
    
		if (cost_alpha-epsilon > (current_loss + c1*alpha*phi_prime_zero)) | (cost_alpha-epsilon >= cost_alphalo)
			alpha_hi = alpha;
			alpha_star = alpha

		else
			grad_line =A*(x + alpha*p)-b;
			
			if abs(grad_line'*p) <= -1*c2*phi_prime_zero-epsilon
				alpha_star = alpha
			break
			end 
 
			if grad_line'*p*(alpha_hi-alpha_lo) >= 0
				alpha_hi = alpha_lo;
				alpha_lo = alpha;
				alpha_star = alpha
			end  
		end  
     	
 	end
   
end   
