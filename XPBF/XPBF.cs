using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Runtime;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using XPBF.SpatialHash;

namespace XPBF
{
    public class XPBF
    {
        //Set up
        private float minX;
        private float minY;
        private float maxX;
        private float maxY;
        float bounce = 0.5f;
        private float g;
        private float dt;
        private float dt2;
        private float gdt2;
        const float PI = (float)Math.PI;

        //Artificial preassure
        private float K;
        private float W_deltaQ;
        private float N; //usually 4-8

        //Constraints count = particles count
        private int count;
        
        //Particles
        private float[] prevX;
        private float[] prevY;
        public float[] posX;
        public float[] posY;
        private SpatialHash2D grid;

        //Constraint
        private float[] invRho0; //rest density
        private float[] alpha;
        private float[] beta;

        private float[] lambda;
        private float[] deltaLambda;

        //Poly6 kernel
        private float poly6H;
        private float poly6H2;
        private float poly6Coeff;

        //Spiky kernel
        private float spikyH;
        private float spikyH2;
        private float spikyCoeff;

        //Query parallelism
        private int threadCount;
        private int chunk;
        private int remainder;
        private List<int>[] perThreadBuffer;

        static int NextPowerOfTwo(int x)
        {
            //bit mask and bit shift to get next smallest power of 2
            x--;
            x |= x >> 1;
            x |= x >> 2;
            x |= x >> 4;
            x |= x >> 8;
            x |= x >> 16;
            x++;
            return x;
        }

        public static (float positionX, float positionY, float restDensity, float reinforcement, float compliance)[]PositionsToParticles(
        double[] x,
        double[] y,
        double restDensity,
        double reinforcement,
        double compliance
        )
        {
            if (x.Length != y.Length)
                throw new ArgumentException("x and y must have same length");

            int count = x.Length;

            float rd = (float)restDensity;
            float d = (float)reinforcement;
            float c = (float)compliance;

            var particles = new (float, float, float, float, float)[count];

            for (int i = 0; i < count; i++)
            {
                particles[i] = (
                    (float)x[i],
                    (float)y[i],
                    rd,
                    d,
                    c
                );
            }

            return particles;
        }

        public XPBF(
            float minX, float minY,
            float maxX, float maxY,
            float gravity,
            float deltaTime,
            float smoothingRadius,

            (float positionX, float positionY, float restDensity, float reinforcement, float compliance)[] particles,

            float artificialPressureK,
            float artificialPressureN,
            float deltaQFactor
        )
        {
            //Domain
            this.minX = minX;
            this.minY = minY;
            this.maxX = maxX;
            this.maxY = maxY;
            this.dt = deltaTime;
            this.dt2 = dt * dt;
            this.g = gravity;
            this.gdt2 = g * dt2;

            //Count
            this.count = particles.Length;

            //Artificial pressure
            this.K = artificialPressureK;
            this.N = artificialPressureN;

            //Allocate particle arrays
            prevX = new float[count];
            prevY = new float[count];
            posX = new float[count];
            posY = new float[count];

            //Constraint arrays
            invRho0 = new float[count];
            alpha = new float[count];
            beta = new float[count];
            lambda = new float[count];
            deltaLambda = new float[count];

            //Init
            for (int i = 0; i < count; i++)
            {
                invRho0[i] = 1f / particles[i].restDensity;
                alpha[i] = particles[i].compliance / dt2;
                beta[i] = particles[i].reinforcement;
                prevX[i] = particles[i].positionX;
                prevY[i] = particles[i].positionY;
                posX[i] = particles[i].positionX;
                posY[i] = particles[i].positionY;
            }

            //Kernel Precomputation
            float h = smoothingRadius;
            float h2 = h * h;
            float h4 = h2 * h2;
            float h6 = h4 * h2;
            float h9 = h6 * h2 * h;

            //Poly6
            poly6H = h;
            poly6H2 = h2;
            poly6Coeff = 315f / (64f * PI * h9);

            //Spiky
            spikyH = h;
            spikyH2 = h2;
            spikyCoeff = -45f / (PI * h6);

            //Artificial pressure setup
            float deltaQ = deltaQFactor * h;
            float deltaQ2 = deltaQ * deltaQ;
            float diff = poly6H2 - deltaQ2;
            W_deltaQ = poly6Coeff * diff * diff * diff;

            //Spatial hash grid
            //bucketCount ≈ 2x–4x active cells count and power of 2(reduces collisions)
            int bucketCount = NextPowerOfTwo(count * 2);
            int capacity = count;
            float cellSize = h;
            grid = new SpatialHash2D(capacity, bucketCount, cellSize);

            //Init threads
            threadCount = 32;
            chunk = count / threadCount;
            remainder = count - chunk * threadCount;
            perThreadBuffer = new List<int>[threadCount];
            for(int m=0;m< threadCount;m++)
            {
                perThreadBuffer[m] = new List<int>(64);
            }
        }

        public void Step(int iterations)
        {
            Array.Clear(lambda, 0, count);
            Predict();

            for (int iter = 0; iter < iterations; iter++)
            {
                ComputeLambda();
                ComputePosCorrection();
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        void Predict()
        {
            grid.Clear();
            for (int i = 0; i < count; i++)
            {
                float px = posX[i];
                float py = posY[i];

                float vx = px - prevX[i];
                float vy = py - prevY[i];

                prevX[i] = px;
                prevY[i] = py;

                px += vx;
                py += vy + gdt2;

                posX[i] = px;
                posY[i] = py;

                grid.Insert(px, py);
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        void ComputeLambda()
        {
            Parallel.For(0, threadCount, t =>
            {
                int start = t * chunk;
                int end = start + chunk;

                if (t == threadCount - 1)
                    end += remainder;

                for (int i = start; i < end; i++)
                {
                    float px = posX[i];
                    float py = posY[i];
                    float invRho = invRho0[i];
                    float a = alpha[i];
                    float b = beta[i];

                    var neighbors = perThreadBuffer[t];
                    grid.Query(px, py, 1, neighbors);
                   
                    //Constraint
                    //sum of poly6
                    float sigmaW = 0;
                    //sum of gradient and its square
                    float sigmaGi_x = 0;
                    float sigmaGi_y = 0;
                    float sigmaGj2 = 0;
                    for (int n = 0; n < neighbors.Count; n++)
                    {
                        int j = neighbors[n];
                        
                        // constraint computation
                        var dx = px - posX[j];
                        var dy = py - posY[j];
                        sigmaW += Poly6(dx, dy);

                        SpikyGrad(dx,dy,out float gx,out float gy);
                        sigmaGi_x += gx;
                        sigmaGi_y += gy;
                        sigmaGj2 += (gx * gx + gy * gy);
                    }
                    float constraint = invRho * sigmaW - 1;
                    float denom = invRho * invRho * (sigmaGj2 + (sigmaGi_x * sigmaGi_x + sigmaGi_y * sigmaGi_y)) + a;
                    float deltalambda = (-constraint - lambda[i] * a * b) / denom;
                    deltaLambda[i] = deltalambda;
                }
            });
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        void ComputePosCorrection()
        {
            Parallel.For(0, threadCount, t =>
            {
                int start = t * chunk;
                int end = start + chunk;

                if (t == threadCount - 1)
                    end += remainder;

                for (int i = start; i < end; i++)
                {
                    //compute deltaPosition
                    float px = posX[i];
                    float py = posY[i];
                    float invRho = invRho0[i];
                    float dl = deltaLambda[i];

                    var neighbors = perThreadBuffer[t];
                    grid.Query(px, py, 1, neighbors);

                    float deltaXx = 0;
                    float deltaXy = 0;
                    for (int n = 0; n < neighbors.Count; n++)
                    {
                        int j = neighbors[n];
                        if (j == i) continue;

                        var dx = px - posX[j];
                        var dy = py - posY[j];

                        //artificial preassure
                        float factor = Poly6(dx, dy)/ W_deltaQ;
                        float scorr = -K;
                        for (int k = 0; k < N; k++)
                        {
                            scorr *= factor;
                        }

                        SpikyGrad(dx, dy, out float gx, out float gy);
                        deltaXx += invRho * gx * (dl + deltaLambda[j] + scorr);
                        deltaXy += invRho * gy * (dl + deltaLambda[j] + scorr);
                    }
                    px += deltaXx;
                    py += deltaXy; // 0 = no bounce, 1 = perfect reflection

                    // X MIN
                    if (px < minX)
                    {
                        px = minX;

                        float vx = px - prevX[i];
                        vx = -vx * bounce;

                        prevX[i] = px - vx;
                    }
                    // X MAX
                    if (px > maxX)
                    {
                        px = maxX;

                        float vx = px - prevX[i];
                        vx = -vx * bounce;

                        prevX[i] = px - vx;
                    }
                    // Y MIN
                    if (py < minY)
                    {
                        py = minY;

                        float vy = py - prevY[i];
                        vy = -vy * bounce;

                        prevY[i] = py - vy;
                    }
                    // Y MAX
                    if (py > maxY)
                    {
                        py = maxY;

                        float vy = py - prevY[i];
                        vy = -vy * bounce;

                        prevY[i] = py - vy;
                    }
                    posX[i] = px;
                    posY[i] = py;
                    lambda[i] += dl;
                }
            });
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        float Poly6(float dx, float dy)
        {
            float r2 = dx * dx + dy * dy;

            if (r2 >= poly6H2) return 0f;

            float diff = poly6H2 - r2;
            return poly6Coeff * diff * diff * diff;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        void SpikyGrad(float dx, float dy, out float gx, out float gy)
        {
            float r2 = dx * dx + dy * dy;

            if (r2 < 1e-12f || r2 >= spikyH2)
            {
                gx = 0f;
                gy = 0f;
                return;
            }

            float invR = 1.0f / (float)Math.Sqrt(r2); // or InvSqrt(r2)
            float r = r2 * invR;

            float diff = spikyH - r;
            float factor = spikyCoeff * diff * diff * invR;

            gx = factor * dx;
            gy = factor * dy;
        }
    }
}
