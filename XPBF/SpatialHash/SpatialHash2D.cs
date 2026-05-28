using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace XPBF.SpatialHash
{
    internal class SpatialHash2D
    {

        private readonly float invCellSize;
        private readonly int bucketCount;/* bucketCount ≈ 2–4×particle count (or active cells) */
        private readonly int capacity;

        private readonly List<int>[] buckets;

        private readonly float[] posX;
        private readonly float[] posY;

        private readonly long[] cellKey;

        private int count;

        public SpatialHash2D(int capacity = 10000, int bucketCount = 16384, float cellSize = 10f)
        {
            if ((bucketCount & (bucketCount - 1)) != 0)
                throw new Exception("bucketCount must be power of 2");
            this.bucketCount = bucketCount;

            this.capacity = capacity;
            this.invCellSize = 1f / cellSize;

            buckets = new List<int>[bucketCount];
            for (int i = 0; i < buckets.Length; i++)
            {
                buckets[i] = new List<int>(4);
            }

            posX = new float[capacity];
            posY = new float[capacity];

            cellKey = new long[capacity];
        }



        public void Insert(float x, float y)
        {
            if (count >= capacity)
                throw new Exception("SpatialHash full");

            //Get Cell Inline
            //use inverse cell size to accelerate computing
            int cx = (int)(x * invCellSize);
            int cy = (int)(y * invCellSize);

            //Hash Function Inline
            int bucketIndex;
            unchecked
            {
                //works better when bucket count is power of 2!!!
                int h = cx * 73856093 ^ cy * 19349663;
                h ^= h >> 13;
                h *= 1013904223;
                h ^= h >> 16;
                bucketIndex = h & (bucketCount - 1);
            }

            int i = count;
            count = i + 1;

            posX[i] = x;
            posY[i] = y;
            cellKey[i] = ((long)cx << 32) | (uint)cy;
            buckets[bucketIndex].Add(i);
        }

        public void Query(float x, float y, int cellRadius, List<int> objectsBuffer)
        {
            objectsBuffer.Clear();

            var bucketsLocal = buckets;
            var keyLocal = cellKey;
            var output = objectsBuffer;

            int cx = (int)(x * invCellSize);
            int cy = (int)(y * invCellSize);

            int minX = cx - cellRadius;
            int maxX = cx + cellRadius;
            int minY = cy - cellRadius;
            int maxY = cy + cellRadius;

            for (int ox = minX; ox <= maxX; ox++)
            {
                //Hash Function Inline
                int hx = ox * 73856093;

                for (int oy = minY; oy <= maxY; oy++)
                {
                    int h;
                    unchecked
                    {
                        h = hx ^ (oy * 19349663);
                        h ^= h >> 13;
                        h *= 1013904223;
                        h ^= h >> 16;
                    }

                    var list = bucketsLocal[h & (bucketCount - 1)];

                    long key = ((long)ox << 32) | (uint)oy;

                    for (int k = 0; k < list.Count; k++)
                    {
                        int i = list[k];
                        if (keyLocal[i] == key)
                        {
                            output.Add(i);
                        }
                    }
                }
            }
        }

        public void Clear()
        {
            for (int i = 0; i < buckets.Length; i++)
            {
                buckets[i].Clear();
            }
            count = 0;
        }
    }
}
