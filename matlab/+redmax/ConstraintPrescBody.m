classdef ConstraintPrescBody < redmax.Constraint
	% ConstraintPrescBody Prescribes body motion.
	
	%%
	properties
		body
		prows % Prescribed motion rows
		vel   % velocity-level integrator?
		q
		qdot
		qddot
	end
	
	%%
	methods
		%%
		function this = ConstraintPrescBody(body,prows,vel)
			this = this@redmax.Constraint(length(prows),0,0,0);
			this.name = [body.name,'-presc'];
			body.presc = this;
			this.body = body;
			this.prows = prows;
			this.vel = vel;
			this.q = zeros(6,1);
			this.qdot = zeros(6,1);
			this.qddot = zeros(6,1);
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [Gm,Gmdot,gm,gmdot,gmddot] = computeJacEqM_(this,Gm,Gmdot,gm,gmdot,gmddot)
			rows = this.idxEM;
			cols = this.body.idxM;
			this.idxQ = cols;
			I = eye(6);
			Gm(rows,cols) = -I(this.prows,:);
			if this.vel
				gmdot(rows) = this.qdot(this.prows);
			else
				gmdot(rows) = this.body.phi(this.prows);
				gmddot(rows) = this.qddot(this.prows);
			end
		end
		
		%%
		function scatterForceEqM_(this,Gmt,lm) %#ok<INUSD>
			this.body.wext_i = this.fcon;
		end
	end
end
