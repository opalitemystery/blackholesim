from flask import Flask, jsonify, request
from flask_cors import CORS
import sys
sys.path.insert(0, 'src')
from runSimulation import runSingleSimulation

app = Flask(__name__)
CORS(app)

def parse_float(value, default=1.0):
    
    try:
        return float(value)
    except (TypeError, ValueError):
        return default

@app.route('/api/simulate', methods=['POST'])
def simulateEndpoint():
    import sys
    
    try:
        print("[BACKEND] Received simulation request", flush=True)
        sys.stdout.flush()
        
        data = request.json
        
        M = parse_float(data.get('M', 1.0))
        a = parse_float(data.get('a', 0.0))
        r0 = parse_float(data.get('r0', 5.0))
        phi0 = parse_float(data.get('phi0', 0.0))
        rDot0 = parse_float(data.get('rDot0', 0.0))
        phiDot0 = parse_float(data.get('phiDot0', 0.192))
        lambdaMax = parse_float(data.get('lambdaMax', 100.0))
        stepInit = parse_float(data.get('stepInit', 0.01))
        outputEvery = int(data.get('outputEvery', 10))
        
        print(f"[BACKEND] Params: M={M}, a={a}, r0={r0}", flush=True)
        sys.stdout.flush()
        
        trajectory, status = runSingleSimulation(
            M, a, r0, phi0, rDot0, phiDot0, lambdaMax, stepInit, outputEvery
        )
        
        print(f"[BACKEND] Got {len(trajectory)} points, status={status}", flush=True)
        sys.stdout.flush()
        
        print("\n" + "="*80, flush=True)
        print(f"TRAJECTORY: M={M}, a={a}, r0={r0}, phiDot0={phiDot0}", flush=True)
        print(f"Status: {status}, Points: {len(trajectory)}", flush=True)
        print("="*80, flush=True)
        for i, point in enumerate(trajectory[:400]): 
            x = point['x']
            y = point['y']
            r = (x**2 + y**2)**0.5
            print(f"  [{i:5d}] x={x:12.8f}, y={y:12.8f}  (r={r:.8f})", flush=True)

        print("="*80 + "\n", flush=True)
        sys.stdout.flush()
        
        return jsonify({
            "success": True,
            "status": status,
            "trajectory": trajectory
        })
    except Exception as e:
        print(f"[ERROR] {str(e)}", flush=True)
        import traceback
        traceback.print_exc()
        return jsonify({"success": False, "error": str(e)}), 400

if __name__ == '__main__':
    print("Starting Kerr Black Hole Geodesics Server...")
    print("Running on http://localhost:5000")
    app.run(debug=True, port=5000)