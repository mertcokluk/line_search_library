%The algorithm is based on Numerical Optimization chapter3 p.56
function [alpha_star] = interpolate(step_len, current_cost, phi_prime_zero, phi_steplen, x, d, ref, order)
	
	if order == 2
		alpha_1 = ((-0.5)*phi_prime_zero*(step_len^2))/(phi_steplen-current_cost-step_len*phi_prime_zero);
		alpha_star=alpha_1;
	end

	if order == 3
		alpha_1 = (-0.5)*phi_prime_zero*(step_len^2)/(phi_steplen-current_cost-phi_prime_zero*step_len);
		z = (1/(((step_len*alpha_1)^2)*(alpha_1-step_len)))*[step_len^2, alpha_1^2; -1*step_len^3, alpha_1^3]*[J((x + alpha_1*d), ref)-current_cost-phi_prime_zero*alpha_1; phi_steplen-current_cost-phi_prime_zero*step_len];
		a = z(1);
		b = z(2);		
		alpha_2 = ( (-1*b)+sqrt(b^2-3*a*phi_prime_zero) )/(3*a);
		alpha_star=alpha_2;
	end

end
