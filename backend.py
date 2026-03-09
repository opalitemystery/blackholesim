from flask import Flask, jsonify, request
from flask_cors import CORS
import sys
sys.path.insert(0, 'src')
from runSimulation import runSingleSimulation, runBatchSimulation

app = Flask(__name__)
CORS(app)

@app.route('/api/simulate', methods=['POST'])
def simulateEndpoint():
    """Single trajectory endpoint"""
    try:
        data = request.json
        
        M = float(data.get('M', 1.0))
        r0 = float(data.get('r0', 5.0))
        phi0 = float(data.get('phi0', 0.0))
        rDot0 = float(data.get('rDot0', 0.0))
        phiDot0 = float(data.get('phiDot0', 0.192))
        lambdaMax = float(data.get('lambdaMax', 100.0))
        stepInit = float(data.get('stepInit', 0.01))
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

@app.route('/api/batch-simulate', methods=['POST'])
def batchSimulateEndpoint():
    """Multiple trajectories endpoint"""
    try:
        data = request.json
        
        M = float(data.get('M', 1.0))
        r0 = float(data.get('r0', 5.0))
        phi0 = float(data.get('phi0', 0.0))
        rDot0 = float(data.get('rDot0', 0.0))
        phiDot0Array = data.get('phiDot0Array', [0.192])
        lambdaMax = float(data.get('lambdaMax', 100.0))
        stepInit = float(data.get('stepInit', 0.01))
        outputEvery = int(data.get('outputEvery', 10))
        
        trajectories = runBatchSimulation(
            M, r0, phi0, rDot0, phiDot0Array, lambdaMax, stepInit, outputEvery
        )
        
        return jsonify({
            "success": True,
            "trajectories": trajectories
        })
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 400

if __name__ == '__main__':
    app.run(debug=True, port=5000)