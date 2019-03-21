%Algorithm is as defined in the book 'Numerical Optimization' p.59
function [step_len] = wolfe(x, grad, A, b, iter)
	
	c1=0.01		%optimization parameters
	c2=0.99;

	p = -1.0*grad;  		 %descent direction
	step_len_old = 0;  %initial previous step length
	step_len = 0.7;      %initial step length 
	step_len_max = 1;  %max step length
	
	epsilon = 0.0001;

	current_cost = 0.5*x'*A*x - b'*x;
	phi_prime_zero = grad'*p;
	phi_steplen = 0.5*(x + step_len*p)'*A*(x + step_len*p) - b'*(x + step_len*p);

	if iter>0	%whether plot step length search process or not
		for s=1:21
			phi(s) = 0.5*(x + ((s-1)/10)*p)'*A*(x + ((s-1)/10)*p) - b'*(x + ((s-1)/10)*p);	 
		end
		mkdir('steplen-analysis');
		fig = plot((0:20)/10, phi, 'b', step_len, phi_steplen, '*');
		title(strcat('wolfe-', 'iteration:', num2str(iter)))
		xlabel('step length')
		ylabel('cost')
		saveas(gcf, strcat('steplen-analysis/', num2str(iter, '%03d'), '-00', '-phi.png'))
	end
	
	trial = 1;
	while step_len-epsilon < step_len_max
		if [ phi_steplen-epsilon > (current_cost + c1*step_len*phi_prime_zero) ] | [ ( phi_steplen >= 0.5*(x + step_len_old*p)'*A*(x + step_len_old*p) - b'*(x + step_len_old*p)) & (step_len_old>0) ]  
			temp = step_len;
			step_len = zoom(x, grad, A, b, c1, c2, step_len_old, step_len); 	%new step length via 'zoom'
			step_len_old = temp;
			phi_steplen = 0.5*(x + step_len*p)'*A*(x + step_len*p) - b'*(x + step_len*p);

			if iter>0	%whether plot step length search process or not
				fig = plot((0:20)/10, phi, 'b', step_len, phi_steplen, '*');
				title(strcat('wolfe-', 'iteration:', num2str(iter)))
				xlabel('step length')
				ylabel('cost')
				saveas(gcf, strcat('steplen-analysis/', num2str(iter, '%03d'), '-', num2str(trial, '%02d'), '-phi.png'))
			end
			trial = trial+1; 
			break
		end

		grad_line = A*(x + step_len*p)-b;
		phi_prime_grad_line = grad'*grad_line;

		if abs(phi_prime_grad_line)  <=  -1*c2*phi_prime_zero 
			break
		end
		if phi_prime_grad_line >= 0
			temp = step_len;
			step_len = zoom(x, grad, A, b, c1, c2, step_len_old, step_len);
			step_len_old = temp;
			phi_steplen = 0.5*(x + step_len*p)'*A*(x + step_len*p) - b'*(x + step_len*p);

			if iter>0	%whether plot step length search process or not
				fig = plot((0:20)/10, phi, 'b', step_len, phi_steplen, '*');
				title(strcat('wolfe-', 'iteration:', num2str(iter)))
				xlabel('step length')
				ylabel('cost')
				saveas(gcf, strcat('steplen-analysis/', num2str(iter, '%03d'), '-', num2str(trial, '%02d'), '-phi.png'))
			end
			trial = trial+1;
			break
		end
 		
		step_len_old = step_len;
		step_len = 0.5*(step_len + step_len_max);
		phi_steplen = 0.5*(x + step_len*p)'*A*(x + step_len*p) - b'*(x + step_len*p);

		if iter>0	%whether plot step length search process or not
				fig = plot((0:20)/10, phi, 'b', step_len, phi_steplen, '*');
				title(strcat('wolfe-', 'iteration:', num2str(iter)))
				xlabel('step length')
				ylabel('cost')
				saveas(gcf, strcat('steplen-analysis/', num2str(iter, '%03d'), '-', num2str(trial, '%02d'), '-phi.png'))
		end
		trial = trial+1; 
	end	
end
