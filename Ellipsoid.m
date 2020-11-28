classdef Ellipsoid < NivoEllipsoid
    properties (Dependent = true, Access = public)
        SemimajorAxis double
        SemiminorAxis double
        InverseFlattening double
        Eccentricity double
    end
    properties (GetAccess = public, SetAccess = private)
        SecondEccentricity = 0;
        Flattening = 0;
        ThirdFlattening = 0;
    end
    properties (Hidden = true, Access = protected)
        a = 1;
        b = 1;
        invf = inf;
        ecc = 0;
        ecc2 = 0;
    end
    properties (GetAccess = protected, Constant = true)
        AngularVelocity = 7292115e-11;
        EarthGravity = 3986004418e+8;
    end

    % GET METHODS
    methods
        function a = get.SemimajorAxis(this)
            a = this.a;
        end
        function b = get.SemiminorAxis(this)
            b = this.b;
        end
        function invf = get.InverseFlattening(this)
            invf = this.invf;
        end
        function ecc = get.Eccentricity(this)
            ecc = this.ecc;
        end
    end
    
    % SET METHODS
    methods
        function this = set.SemimajorAxis(this, a)
            this.a = a;
            this.b = (1 - this.Flattening) * a;
        end
        function this = set.SemiminorAxis(this, b)
            a_ = this.a;
            this.b = b;
            this.ecc = sqrt(a_*a_ - b*b) / a_;
            this.ecc2 = sqrt(a_*a_ - b*b) / b;
            f = (a_ - b) / a_;
            this.invf = 1 / f;
            this.Flattening = f;
            this.ThirdFlattening = f / (2 - f);
            this.SecondEccentricity = this.ecc2;
        end
        function this = set.InverseFlattening(this, invf)
            f = 1 / invf;
            this.ecc = sqrt((2 - f) * f);
            this.ecc2 = this.ecc / (1 - f);
            this.invf = invf;
            this.Flattening = f;
            this.ThirdFlattening = f / (2 - f);
            this.b = (1 - this.Flattening) * this.a;
            this.SecondEccentricity = this.ecc2;
        end
        function this = set.Eccentricity(this, ecc)
            this.ecc = ecc;
            e2 = ecc * ecc;
            this.b = (1 - f) * this.a;
            this.ecc2 = e2 - (1 - e2);
            this.invf = 1 / f;
            this.Flattening = f;
            this.ThirdFlattening = f / (2 - f);
            this.SecondEccentricity = this.ecc2;
        end
    end
    
    methods
        %   args:
        %
        %
        %   returns:
        %
        %
        function [x, y, z] = geog2geoc(ellipsoid, B, L, latType, h)
            if nargin < 4
                latType = 'geographic';
            end
            if nargin < 5; h = 0; end
            
            a_ = ellipsoid.SemimajorAxis;
            b_ = ellipsoid.SemiminorAxis;
            e = ellipsoid.Eccentricity; e2 = e * e;
            eu = ellipsoid.SecondEccentricity; eu2 = eu * eu;
            
            sinf = @sind;
            cosf = @cosd;
            tanf = @tand;
            atanf = @atand;
            
            switch latType
                case 'reduced'
                    B = atanf((a_/b_) * tanf(B));
                case 'geocentric'
                    B = atanf((a_*a_ / b_*b_) * tanf(B));
            end
            
            n2 = eu2 * cosf(B).^2;
            V = sqrt(1 + n2);
            W = sqrt(1 - e2 * sinf(B).^2);
            
            xv = (a_*a_ ./ (b_*V)+h) .* cosf(B) .* cosf(L);
            yv = (a_*a_ ./ (b_*V)+h) .* cosf(B) .* sinf(L);
            zv = (b_./V+h) .* sinf(B);
            
            xw = (a_./W+h) .* cosf(B) .* cosf(L);
            yw = (a_./W+h) .* cosf(B) .* sinf(L);
            zw = (b_*b_ ./ (a_*W)+h) .* sinf(B);
            
            x = round((xv + xw) / 2, 3);
            y = round((yv + yw) / 2, 3);
            z = round((zv + zw) / 2, 3);
        end
        
        %   args:
        %
        %
        %   returns:
        %
        %
        function [b, l, h] = geoc2geog(ellipsoid, x, y, z)
            a_ = ellipsoid.SemimajorAxis;
            b_ = ellipsoid.SemiminorAxis;
            eu = ellipsoid.SecondEccentricity; eu2 = eu * eu;
            
            atanf = @atand;
            atan2f = @atan2d;
            
            x_ = x; y_ = y; z_ = z;
            z__ = (b_*z_) / sqrt(x_*x_ + y_*y_ + z_*z_);
            k_ = (a_ / b_) * (a_ / b_);
            z_old = inf;
            while abs(z__ - z_old) > 1e-7
                z_old = z__;
                k = sqrt(x_*x_ + y_*y_ + (z_ + eu2*z__)*(z_ + eu2*z__));
                l_ = x_ / k; m = y_ / k; n = (z_ + eu2*z__) / k;
                r = l_*l_ + m*m + k_*n*n;
                s = l_*x_ + m*y_ + k_*n*z_;
                t = x_*x_ + y_*y_ + k_*z_*z_ - a_*a_;
                h = (s + sqrt(s*s - r*t)) / r;
                root = sqrt(s*s - r*t);
                
                z__ = z_ - n*((s - root) / r);
                if h < 0; z__ = z_ - n*h; end
            end
            h = (s - root) / r;
            if h < 0; h = (s + root) / r; end
            
            b = atanf((z_ + eu2*z__) / sqrt(x_*x_ + y_*y_));
            l = atan2f(y_, x_);
        end
        
        %   args:
        %
        %
        %   returns:
        %
        %
        function [A1, A2, S] = JTP2(ellipsoid, B1, L1, B2, L2)
            b_ = ellipsoid.SemiminorAxis;
            f = ellipsoid.Flattening;
            ec = ellipsoid.SecondEccentricity;
            ec2 = ec * ec;
            
            B1 = deg2rad(B1); B2 = deg2rad(B2);
            L1 = deg2rad(L1); L2 = deg2rad(L2); 
            
            u = atan((1 - f) * tan(B1));
            u_ = atan((1 - f) * tan(B2));
            w = L2 - L1 ; dL = w;
            
            dL_old = inf;
            while abs(dL - dL_old) > 1e-12
                dL_old = dL;
                snq = (sqrt((cos(u_) * sin(dL))^2 + (cos(u) * sin(u_) - sin(u) * cos(u_) * cos(dL))^2));
                csq = sin(u) * sin(u_) + cos(u) * cos(u_) * cos(dL);
                
                q = atan2(snq, csq);
                
                Al = asin(cos(u) * cos(u_) * sin(dL) / snq);
                km = csq - (2 * sin(u) * sin(u_) / cos(Al)^2);
                c = (f / 16) * cos(Al)^2 * (4 + f * (4 - 3 * cos(Al)^2));
                dL = w + (1 - c) * f * sin(Al) * (q + c * snq * (km + c * csq * (-1 + 2 * km^2)));
            end
            
            u2_ = ec2 * cos(Al)^2;
            a__ = 1 + (u2_ / 16384) * (4096 + u2_ * (-768 + u2_ * (320 - 175 * u2_)));
            b__ = (u2_ / 1024) * (256 + u2_ * (-128 + u2_ * (74 - 47*u2_)));
            dq = b__ * snq * (km + (b__ / 4) * (csq * (-1 + 2 * km^2) - (b__ / 6) * km * (-3 + 4 * snq^2) * (-3 + 4 * km^2)));
            
            S = b_ * a__ *abs(q - dq);
            A1 = atan2(cos(u_) * sin(dL), (cos(u) * sin(u_) - sin(u) * cos(u_) * cos(dL)));
            A2 = atan2(cos(u) * sin(dL), (-sin(u) * cos(u_) + cos(u) * sin(u_) * cos(dL)));
            if A1 < pi
                A2 = rad2deg(A2 + pi);
            else
                A2 = rad2deg(A2 - pi);
            end
            A1 = rad2deg(A1);
        end
    end
end