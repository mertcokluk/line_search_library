%Algorithm is as defined in the lecture notes BMS Basic Course: Nonlinear Optimization p.17
function [step_len] = armijo(x, grad, A, b, iter)

	nfft_real_sqd = size(x, 3)^2;
	sigma=0.1/nfft_real_sqd; 		%optimization parameters
	order=2;		        %order of the interpolation polynomial

	epsilon = 0.1;
		
	d = -1.0*grad;   		%descent direction
	step_len = 9; 		%initial step length

	current_cost = 0.5*x'*A*x - b'*x;
	phi_prime_zero = sum(sum(sum(conj(grad).*d, 1), 2), 3);
	phi_steplen = 0.5*(x + step_len*d)'*A*(x + step_len*d) - b'*(x + step_len*d);	
	
	if iter>0	%whether plot step length search process or not
		for s=1:111
			phi(s) = J((x + ((s-1)/10)*d), ref);
		end
		mkdir('steplen-analysis');
		fig = plot((0:110)/10, phi, 'b', step_len, phi_steplen, '*');
		title(strcat('armijo-', 'iteration:', num2str(iter)))
		xlabel('step length')
		ylabel('costIVA')
		saveas(gcf, strcat('steplen-analysis/', num2str(iter, '%03d'), '-00', '-phi.png'))
	end
	
	trial = 1;
    	while ( phi_steplen-epsilon > (current_cost + sigma*step_len*phi_prime_zero) ) 			  	%Armijo condition 
		est = interpolate(step_len, current_cost, phi_prime_zero, phi_steplen, x, d, ref, order);	%interpolate for update step length
		
		if est > 0.9*step_len		%Barrier for extremely different interpolation estimate-based step length
			step_len = 0.9*step_len;
		elseif est < 0.7*step_len
			step_len = 0.7*step_len;
		else
			step_len = est;
		end

		phi_steplen = 0.5*(x + step_len*p)'*A*(x + step_len*p) - b'*(x + step_len*p);
		c = phi_steplen-epsilon - (current_cost + sigma*step_len*phi_prime_zero)
		
		if iter>0 			%whether plot step length search process or not
			fig = plot((0:110)/10, phi, 'b', step_len, phi_steplen, '*');
			title(strcat('armijo-', 'iteration:', num2str(iter)))
			xlabel('step length')
			ylabel('costIVA')
			saveas(gcf, strcat('steplen-analysis/', num2str(iter, '%03d'), '-', num2str(trial, '%02d'), '-phi.png'))
		end
		trial = trial+1; 
    	end

end
