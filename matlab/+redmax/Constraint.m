classdef (Abstract) Constraint < handle
	% Constraint Generic constraint
	% A constraint can be applied in reduced or maximal coordinates.
	% Reduced: keeping a quaternion to be of unit length.
	% Maximal: Holding two bodies together in world space.
	
	%%
	properties
		name
		nconEM  % Number of maximal equality constraints
		nconER  % Number of reduced equality constraints
		nconIM  % Number of maximal inequality constraints
		nconIR  % Number of reduced inequality constraints
		idxEM   % Maximal equality constraint indices
		idxER   % Reduced equality constraint indices
		idxIM   % Maximal inequality constraint indices
		idxIR   % Reduced inequality constraint indices
		idxQ    % Associated DOF indices
		idxT    % Tangent matrix indices
		activeM % Whether the maximal inequality constraint is active
		activeR % Whether the reduced inequality constraint is active
		fcon    % Computed constraint force
		contacts% Contact points for loop-closing constraints
		next    % Next constraint in traversal order
	end
	
	%%
	methods
		%%
		function this = Constraint(nconEM,nconER,nconIM,nconIR)
			this.nconEM = nconEM;
			this.nconER = nconER;
			this.nconIM = nconIM;
			this.nconIR = nconIR;
			this.activeM = false;
			this.activeR = false;
			this.contacts = {};
		end
		
		%%
		function countDofs(this)
			% Counts DOFs
			nem = redmax.Scene.countEM();
			ner = redmax.Scene.countER();
			nim = redmax.Scene.countIM();
			nir = redmax.Scene.countIR();
			this.idxEM = nem + (1:this.nconEM);
			this.idxER = ner + (1:this.nconER);
			this.idxIM = nim + (1:this.nconIM);
			this.idxIR = nir + (1:this.nconIR);
			nem = nem + this.nconEM;
			ner = ner + this.nconER;
			nim = nim + this.nconIM;
			nir = nir + this.nconIR;
			redmax.Scene.countEM(nem);
			redmax.Scene.countER(ner);
			redmax.Scene.countIM(nim);
			redmax.Scene.countIR(nir);
		end
		
		%%
		function [listM,listR] = getActiveList(this,listM,listR)
			% Gets list of active inequality indices
			if nargin == 1
				listM = [];
				listR = [];
			end
			if this.activeM
				listM(end+1) = this.idxIM;
			end
			if this.activeR
				listR(end+1) = this.idxIR;
			end
			if ~isempty(this.next)
				[listM,listR] = this.next.getActiveList(listM,listR);
			end
		end
		
		%%
		function [Gm,Gmdot,gm,gmdot,gmddot] = computeJacEqM(this,Gm,Gmdot,gm,gmdot,gmddot)
			if nargin == 1
				nm = redmax.Scene.countM();
				nem = redmax.Scene.countEM();
				Gm = zeros(nem,nm);
				Gmdot = zeros(nem,nm);
				gm = zeros(nem,1);
				gmdot = zeros(nem,1);
				gmddot = zeros(nem,1);
			end
			[Gm,Gmdot,gm,gmdot,gmddot] = this.computeJacEqM_(Gm,Gmdot,gm,gmdot,gmddot);
			if ~isempty(this.next)
				[Gm,Gmdot,gm,gmdot,gmddot] = this.next.computeJacEqM(Gm,Gmdot,gm,gmdot,gmddot);
			end
		end
		
		%%
		function [Gr,Grdot,gr,grdot,grddot] = computeJacEqR(this,Gr,Grdot,gr,grdot,grddot)
			if nargin == 1
				nr = redmax.Scene.countR();
				ner = redmax.Scene.countER();
				Gr = zeros(ner,nr);
				Grdot = zeros(ner,nr);
				gr = zeros(ner,1);
				grdot = zeros(ner,1);
				grddot = zeros(ner,1);
			end
			[Gr,Grdot,gr,grdot,grddot] = this.computeJacEqR_(Gr,Grdot,gr,grdot,grddot);
			if ~isempty(this.next)
				[Gr,Grdot,gr,grdot,grddot] = this.next.computeJacEqR(Gr,Grdot,gr,grdot,grddot);
			end
		end
		
		%%
		function [Cm,Cmdot,cm,cmdot,cmddot] = computeJacIneqM(this,Cm,Cmdot,cm,cmdot,cmddot)
			if nargin == 1
				nm = redmax.Scene.countM();
				nim = redmax.Scene.countIM();
				Cm = zeros(nim,nm);
				Cmdot = zeros(nim,nm);
				cm = zeros(nim,1);
				cmdot = zeros(nim,1);
				cmddot = zeros(nim,1);
			end
			[Cm,Cmdot,cm,cmdot,cmddot] = this.computeJacIneqM_(Cm,Cmdot,cm,cmdot,cmddot);
			if ~isempty(this.next)
				[Cm,Cmdot,cm,cmdot,cmddot] = this.next.computeJacIneqM(Cm,Cmdot,cm,cmdot,cmddot);
			end
		end
		
		%%
		function [Cr,Crdot,cr,crdot,crddot] = computeJacIneqR(this,Cr,Crdot,cr,crdot,crddot)
			if nargin == 1
				nr = redmax.Scene.countR();
				nir = redmax.Scene.countIR();
				Cr = zeros(nir,nr);
				Crdot = zeros(nir,nr);
				cr = zeros(nir,1);
				crdot = zeros(nir,1);
				crddot = zeros(nir,1);
			end
			[Cr,Crdot,cr,crdot,crddot] = this.computeJacIneqR_(Cr,Crdot,cr,crdot,crddot);
			if ~isempty(this.next)
				[Cr,Crdot,cr,crdot,crddot] = this.next.computeJacIneqR(Cr,Crdot,cr,crdot,crddot);
			end
		end
		
		%%
		function scatterForceEqM(this,Gmt,lm)
			if this.nconEM > 0
				this.fcon = -Gmt(this.idxQ,this.idxEM)*lm(this.idxEM,1);
			else
				this.fcon = zeros(length(this.idxQ),1);
			end
			this.scatterForceEqM_(Gmt,lm);
			if ~isempty(this.next)
				this.next.scatterForceEqM(Gmt,lm);
			end
		end
		
		%%
		function scatterForceEqR(this,Grt,lr)
			if this.nconER > 0
				this.fcon = -Grt(this.idxQ,this.idxER)*lr(this.idxER,1);
			else
				this.fcon = zeros(length(this.idxQ),1);
			end
			this.scatterForceEqR_(Grt,lr);
			if ~isempty(this.next)
				this.next.scatterForceEqR(Grt,lr);
			end
		end
		
		%%
		function scatterForceIneqM(this,Cmt,lm)
			if this.nconIM > 0
				this.fcon = -Cmt(this.idxQ,this.idxIM)*lm(this.idxIM,1);
			else
				this.fcon = zeros(length(this.idxQ),1);
			end
			this.scatterForceIneqM_(Cmt,lm);
			if ~isempty(this.next)
				this.next.scatterForceIneqM(Cmt,lm);
			end
		end
		
		%%
		function scatterForceIneqR(this,Crt,lr)
			if this.nconIR > 0
				this.fcon = -Crt(this.idxQ,this.idxIR)*lr(this.idxIR,1);
			else
				this.fcon = zeros(length(this.idxQ),1);
			end
			this.scatterForceIneqR_(Crt,lr);
			if ~isempty(this.next)
				this.next.scatterForceIneqR(Crt,lr);
			end
		end
		
		%%
		function [value,isterminal,direction] = ineqEventFcn(this)
			if nargin == 1
				value = [];
				isterminal = [];
				direction = [];
			end
			[value,isterminal,direction] = this.ineqEventFcn_(value,isterminal,direction);
			if ~isempty(this.next)
				[value,isterminal,direction] = this.next.ineqEventFcn_(value,isterminal,direction);
			end
		end
		
		%%
		function ineqProjPos(this)
			this.ineqProjPos_();
			if ~isempty(this.next)
				this.next.ineqProjPos();
			end
		end
		
		%%
		function generateContactsJoint(this)
			% Generates the local contact structures for joint contacts
			this.generateContactsJoint_();
			for i = 1 : length(this.contacts)
				this.contacts{i}.a = 0;
			end
			if ~isempty(this.next)
				this.next.generateContactsJoint();
			end
		end
		
		%%
		function generateContactsCollision(this)
			% Generates the local contact structures for collision contacts
			this.generateContactsCollision_();
			for i = 1 : length(this.contacts)
				this.contacts{i}.a = 0;
			end
			if ~isempty(this.next)
				this.next.generateContactsCollision();
			end
		end
		
		%%
		function computeContactMultiplier(this,h,SPreg)
			% Computes the contact Lagrange multiplier
			this.computeContactMultiplier_(h,SPreg);
			if ~isempty(this.next)
				this.next.computeContactMultiplier(h,SPreg);
			end
		end
		
		%%
		function T = computeTangentMatrix(this,T)
			if nargin == 1
				nt = redmax.Scene.countT();
				nm = redmax.Scene.countM();
				T = zeros(nt,nm);
			end
			T = this.computeTangentMatrix_(T);
			if ~isempty(this.next)
				T = this.next.computeTangentMatrix(T);
			end
		end
		
		%%
		function [bl,bu,idx] = computeFrictionLimits(this,mu,SPathresh,bl,bu,idx)
			if nargin == 3
				nt = redmax.Scene.countT();
				bl = zeros(nt,1);
				bu = zeros(nt,1);
				idx = [];
			end
			[bl,bu,idx] = this.computeFrictionLimits_(mu,SPathresh,bl,bu,idx);
			if ~isempty(this.next)
				[bl,bu,idx] = this.next.computeFrictionLimits(mu,SPathresh,bl,bu,idx);
			end
		end
		
		%%
		function draw(this)
			this.draw_();
			if ~isempty(this.next)
				this.next.draw();
			end
		end
	end
	
	%%
	methods (Access = protected)
		%%
		function [Gm,Gmdot,gm,gmdot,gmddot] = computeJacEqM_(this,Gm,Gmdot,gm,gmdot,gmddot) %#ok<INUSL>
		end
		
		%%
		function [Gr,Grdot,gr,grdot,grddot] = computeJacEqR_(this,Gr,Grdot,gr,grdot,grddot) %#ok<INUSL>
		end
		
		%%
		function [Cm,Cmdot,cm,cmdot,cmddot] = computeJacIneqM_(this,Cm,Cmdot,cm,cmdot,cmddot) %#ok<INUSL>
		end
		
		%%
		function [Cr,Crdot,cr,crdot,crddot] = computeJacIneqR_(this,Cr,Crdot,cr,crdot,crddot) %#ok<INUSL>
		end
		
		%%
		function scatterForceEqM_(this,Gmt,lm) %#ok<INUSD>
		end
		
		%%
		function scatterForceEqR_(this,Grt,lr) %#ok<INUSD>
		end
		
		%%
		function scatterForceIneqM_(this,Cmt,lm) %#ok<INUSD>
		end
		
		%%
		function scatterForceIneqR_(this,Crt,lr) %#ok<INUSD>
		end
		
		%%
		function [value,isterminal,direction] = ineqEventFcn_(this,value,isterminal,direction) %#ok<INUSL>
		end
		
		%%
		function ineqProjPos_(this) %#ok<MANU>
		end
		
		%%
		function generateContactsJoint_(this) %#ok<MANU>
		end
		
		%%
		function generateContactsCollision_(this) %#ok<MANU>
		end
		
		%%
		function computeContactMultiplier_(this,h,SPreg) %#ok<INUSD>
		end
		
		%%
		function T = computeTangentMatrix_(this,T) %#ok<INUSL>
		end
		
		%%
		function [bl,bu,idx] = computeFrictionLimits_(this,mu,SPathresh,bl,bu,idx) %#ok<INUSL>
		end
		
		%%
		function draw_(this) %#ok<MANU>
		end
	end
end
