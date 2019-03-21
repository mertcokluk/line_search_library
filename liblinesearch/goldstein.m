%Algorithm is as defined in the script page about Goldstein-Armijo search
function [step_len] = goldstein(ref, x, grad, iter)

	nfft_real_sqd = size(x, 3)^2;
	mu1=0.1/nfft_real_sqd;	%optimization parameters
	mu2=500/nfft_real_sqd;
	rho1=0.9;
	rho2=1.05;

	epsilon = 0.1;
	
	d = -1.0*grad;   %descent direction
	step_len = 9;  %initial step length
	
	current_cost = J(x, ref);
	phi_prime_zero = sum(sum(sum(conj(grad).*d, 1), 2), 3); 
	phi_steplen = J((x + step_len*d), ref);	
	
	if iter>0	%whether plot step length search process or not
		for s=1:21
			phi(s) = J((x + ((s-1)/10)*d), ref);
		end
		mkdir('steplen-analysis');
		fig = plot((0:20)/10, phi, 'b', step_len, phi_steplen, '*');
		title(strcat('goldstein-', 'iteration:', num2str(iter)))
		xlabel('step length')
		ylabel('costIVA')
		saveas(gcf, strcat('steplen-analysis/', num2str(iter, '%03d'), '-00', '-phi.png'))
	end
	
	trial = 1;	
	while 1>0
		if ( (current_cost + mu2*step_len*phi_prime_zero)-epsilon > phi_steplen )  %condition G.2 is violated
			c = (current_cost + mu2*step_len*phi_prime_zero)-epsilon - phi_steplen		
			step_len = step_len * rho2;	%update step length
			phi_steplen = J((x + step_len*d), ref);

			if iter>0	%whether plot step length search process or not
				fig = plot((0:20)/10, phi, 'b', step_len, phi_steplen, '*');
				title(strcat('goldstein-', 'iteration:', num2str(iter)))
				xlabel('step length')
				ylabel('costIVA')
				saveas(gcf, strcat('steplen-analysis/', num2str(iter, '%03d'), '-', num2str(trial, '%02d'), '-phi.png'))
			end
			trial = trial+1;

		elseif phi_steplen-epsilon > (current_cost + mu1*step_len*phi_prime_zero)	%sufficient decrease, i.e. G.1, is violated
			b = phi_steplen-epsilon - (current_cost + mu1*step_len*phi_prime_zero)
			step_len = step_len * rho1;
			phi_steplen = J((x + step_len*d), ref);

			if iter>0	%whether plot step length search process or not
				fig = plot((0:20)/10, phi, 'b', step_len, phi_steplen, '*');
				title(strcat('goldstein-', 'iteration:', num2str(iter)))
				xlabel('step length')
				ylabel('costIVA')
				saveas(gcf, strcat('steplen-analysis/', num2str(iter, '%03d'), '-', num2str(trial, '%02d'), '-phi.png'))
			end
			trial = trial+1;

		else 			%Both sides of Goldstein rule are fulfilled.
			break
		end	
	end

end
