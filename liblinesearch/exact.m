%Assign the best possible step length  
function [step_len] = exact(x, grad, A, b, iter)
	
	p = -1.0*grad;   %descent direction

	if iter>0	%whether plot step length search process or not

		for s=1:111	%the slice along the descent direction
			phi(s) = 0.5*(x + ((s-1)/50)*p)'*A*(x + ((s-1)/50)*p) - b'*(x + ((s-1)/50)*p); 
		end

		sstar=1;	%argmin of the vector phi
		for t=2:111
			if phi(t)<phi(sstar)
				sstar=t;
			end
		end

		step_len = (sstar-1)/50;
		mkdir('steplen-analysis');
		fig = plot((0:110)/50, phi, 'b', step_len, phi(sstar), '*');
		title(strcat('exact line search-', 'iteration:', num2str(iter)))
		xlabel('step length')
		ylabel('cost')
		saveas(gcf, strcat('steplen-analysis/', num2str(iter, '%03d'), '-00', '-phi.png'))
	end
end
