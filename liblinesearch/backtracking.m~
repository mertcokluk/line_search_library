%Algorithm is as defined in the book 'Numerical Optimization' p.41
function [step_len] = backtracking(x, grad, A, b, iter)
	
	nfft_real_sqd = size(x, 3)^2;
	c= 0.1/nfft_real_sqd;		%optimization paramaters 
	rho=0.9;
	
	p = -1.0*grad;   %descent direction
	step_len = 2;  %initial step length
	
	current_cost = 0.5*x'*A*x - b'*x;
	phi_prime_zero = sum(sum(sum(conj(grad).*p, 1), 2), 3);
	phi_steplen = 0.5*(x + step_len*p)'*A*(x + step_len*p) - b'*(x + step_len*p);
	epsilon = 0.001;

	if iter>0	%whether plot step length search process or not
		for s=1:21
			phi(s) =  0.5*(x + ((s-1)/10)*p)'*A*(x + ((s-1)/10)*p) - b'*(x + ((s-1)/10)*p);	
		end
		mkdir('steplen-analysis');
		fig = plot((0:20)/10, phi, 'b', step_len, phi_steplen, '*');
		title(strcat('backtracking-', 'iteration:', num2str(iter)))
		xlabel('step length')
		ylabel('cost')
		saveas(gcf, strcat('steplen-analysis/', num2str(iter, '%03d'), '-00', '-phi.png'))
	end
	
	trial = 1;
	while phi_steplen-epsilon > ( current_cost + c*step_len*phi_prime_zero )   %sufficient decrease
		d = phi_steplen-epsilon - ( current_cost + c*step_len*phi_prime_zero )
		step_len = step_len * rho;		%update step length
		phi_steplen = 0.5*(x + step_len*p)'*A*(x + step_len*p) - b'*(x + step_len*p);

		if iter>0	%whether plot step length search process or not
			fig = plot((0:20)/10, phi, 'b', step_len, phi_steplen, '*');
			title(strcat('backtracking-', 'iteration:', num2str(iter)))
			xlabel('step length')
			ylabel('cost')
			saveas(gcf, strcat('steplen-analysis/', num2str(iter, '%03d'), '-', num2str(trial, '%02d'), '-phi.png'))
		end
		trial = trial+1; 
	end

end
