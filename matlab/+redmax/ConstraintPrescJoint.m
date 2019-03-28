classdef ConstraintPrescJoint < redmax.Constraint
	% ConstraintPrescJoint Prescribes joint motion.
	
	%%
	properties
		joint
		vel   % velocity-level integrator?
		q
		qdot
		qddot
	end
	
	%%
	methods
		%%
		function this = ConstraintPrescJoint(joint,vel)
			this = this@redmax.Constraint(0,joint.ndof,0,0);
			this.name = [joint.name,'-presc'];
			joint.presc = this;
			this.joint = joint;
			this.vel = vel;
			this.q = zeros(joint.ndof,1);
			this.qdot = zeros(joint.ndof,1);
			this.qddot = zeros(joint.ndof,1);
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [Gr,Grdot,gr,grdot,grddot] = computeJacEqR_(this,Gr,Grdot,gr,grdot,grddot)
			rows = this.idxER;
			cols = this.joint.idxR;
			this.idxQ = cols;
			Gr(rows,cols) = -eye(this.joint.ndof);
			gr(rows) = this.q - this.joint.q;
			if this.vel
				grdot(rows) = this.qdot;
			else
				grdot(rows) = this.joint.qdot;
				grddot(rows) = this.qddot;
			end
		end
		
		%%
		function scatterForceEqR_(this,Grt,lr) %#ok<INUSD>
			this.joint.tau = this.fcon;
		end
	end
end
