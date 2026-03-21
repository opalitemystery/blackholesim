from flask import Flask, jsonify, request
from flask_cors import CORS
import sys
sys.path.insert(0, 'src')
from runSimulation import runSingleSimulation

app = Flask(__name__)
CORS(app)

def parse_float(value, default=1.0):
    """Parse float with high precision, handling string inputs"""
    try:
        return float(value)
    except (TypeError, ValueError):
        return default

@app.route('/api/simulate', methods=['POST'])
def simulateEndpoint():
    """Single trajectory endpoint - high precision"""
    try:
        data = request.json
        
        M = parse_float(data.get('M', 1.0))
        r0 = parse_float(data.get('r0', 5.0))
        phi0 = parse_float(data.get('phi0', 0.0))
        rDot0 = parse_float(data.get('rDot0', 0.0))
        phiDot0 = parse_float(data.get('phiDot0', 0.192))
        lambdaMax = parse_float(data.get('lambdaMax', 100.0))
        stepInit = parse_float(data.get('stepInit', 0.01))
        outputEvery = int(data.get('outputEvery', 10))
        
        trajectory, status = runSingleSimulation(
            M, r0, phi0, rDot0, phiDot0, lambdaMax, stepInit, outputEvery
        )
        
        return jsonify({
            "success": True,
            "status": status,
            "trajectory": trajectory
        })
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 400

if __name__ == '__main__':
    print("Starting Black Hole Geodesics Server...")
    print("Running on http://localhost:5000")
    app.run(debug=True, port=5000)