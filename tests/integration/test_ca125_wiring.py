import pytest
import asyncio
from unittest.mock import MagicMock, AsyncMock
from api.services.resistance_prophet_service import ResistanceProphetService
from api.services.resistance_prophet.schemas import ResistanceSignalData, ResistanceSignal

@pytest.mark.asyncio
async def test_ca125_wiring_absent():
    """Verify that without CA125 Service, the signal is skipped (not evaluated)."""
    service = ResistanceProphetService(ca125_service=None)
    
    # Predict
    result = await service.predict_resistance(
        patient_id="TEST_PATIENT",
        disease_type="ovarian",
        current_sae_features={}, 
        return_object=True
    )
    
    # Assert
    # Signal list should NOT contain CA125_KINETICS
    signal_types = [s.signal_type for s in result.signals_detected]
    assert ResistanceSignal.CA125_KINETICS not in signal_types, "CA125 Signal should not be evaluated when service is None"

@pytest.mark.asyncio
async def test_ca125_wiring_present():
    """Verify that WITH CA125 Service, the signal is evaluated."""
    # Mock the Intelligence Service
    mock_ca125 = MagicMock()
    mock_result = ResistanceSignalData(
        signal_type=ResistanceSignal.CA125_KINETICS,
        detected=True,
        probability=0.99,
        confidence=1.0,
        rationale="Mocked Detected",
        provenance={"status": "mock"},
        baseline_reliability=1.0
    )
    # Async mock for analyze_kinetics
    mock_ca125.analyze_kinetics = AsyncMock(return_value=mock_result)
    
    service = ResistanceProphetService(ca125_service=mock_ca125)
    
    # Predict with Injection
    history = [{"date": "2023-01-01", "value": 10}, {"date": "2023-02-01", "value": 100}]
    result = await service.predict_resistance(
        patient_id="TEST_PATIENT",
        disease_type="ovarian",
        current_sae_features={}, 
        ca125_history_injection=history,
        return_object=True
    )
    
    # Assert
    signal_types = [s.signal_type for s in result.signals_detected]
    assert ResistanceSignal.CA125_KINETICS in signal_types, "CA125 Signal MUST be evaluated when service is injected"
    
    # Find the signal
    sig = next(s for s in result.signals_detected if s.signal_type == ResistanceSignal.CA125_KINETICS)
    assert sig.detected is True
    assert sig.probability == 0.99
    assert sig.rationale == "Mocked Detected"
