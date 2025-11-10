# MODULE 8: BACKEND REFRESH STATUS ENDPOINT

## **Purpose**

Expose `/api/trials/refresh_status` endpoint using existing trial refresh service.

---

## **File Location**

Option A: Add to `oncology-coPilot/oncology-backend/main.py`  
Option B: Create `oncology-coPilot/oncology-backend/backend/routers/trials.py`

---

## **Endpoint Specification**

### **Endpoint:**
```
POST /api/trials/refresh_status
```

### **Request Body:**
```typescript
{
  nct_ids: string[];
  state_filter?: string | null;  // Optional: "NY", "CA", etc.
}
```

### **Response:**
```typescript
{
  refreshed_count: number;
  trial_data: {
    [nctId: string]: {
      status: string;
      locations: Location[];
      last_updated: string;
    };
  };
}
```

---

## **Implementation**

### **Option A: Add to main.py**

```python
from backend.services.trial_refresh import refresh_trial_status_with_retry

@app.post("/api/trials/refresh_status")
async def refresh_status_endpoint(request: Dict[str, Any] = Body(...)):
    """
    Refresh live trial status and locations from ClinicalTrials.gov API v2.
    
    Body: {
        "nct_ids": ["NCT12345", "NCT67890"],
        "state_filter": "NY"  # Optional
    }
    
    Returns:
        {
            "refreshed_count": 2,
            "trial_data": {
                "NCT12345": {
                    "status": "RECRUITING",
                    "locations": [...],
                    "last_updated": "2024-10-20T12:00:00Z"
                }
            }
        }
    """
    nct_ids = request.get("nct_ids", [])
    state_filter = request.get("state_filter")
    
    if not nct_ids:
        raise HTTPException(status_code=400, detail="nct_ids required")
    
    if not isinstance(nct_ids, list):
        raise HTTPException(status_code=400, detail="nct_ids must be a list")
    
    # Refresh trials using existing service
    refreshed_data = await refresh_trial_status_with_retry(nct_ids)
    
    # Apply state filter if requested
    if state_filter:
        from backend.services.trial_refresh.filters import filter_locations_by_state
        filtered_data = {}
        for nct_id, data in refreshed_data.items():
            filtered_locations = filter_locations_by_state(
                data.get("locations", []),
                state_filter
            )
            if filtered_locations:  # Only include if has matching locations
                filtered_data[nct_id] = {
                    **data,
                    "locations": filtered_locations
                }
        refreshed_data = filtered_data
    
    return {
        "refreshed_count": len(refreshed_data),
        "trial_data": refreshed_data
    }
```

---

## **Testing**

### **Manual Test:**
```bash
curl -X POST http://localhost:8000/api/trials/refresh_status \
  -H "Content-Type: application/json" \
  -d '{"nct_ids": ["NCT12345"], "state_filter": "NY"}'
```

### **Expected Response:**
```json
{
  "refreshed_count": 1,
  "trial_data": {
    "NCT12345": {
      "status": "RECRUITING",
      "locations": [
        {
          "facility": "Memorial Sloan Kettering",
          "city": "New York",
          "state": "NY",
          "status": "recruiting",
          "contact_phone": "212-639-XXXX"
        }
      ],
      "last_updated": "2024-10-20T12:00:00Z"
    }
  }
}
```

---

## **Acceptance Criteria**

- [ ] Endpoint accepts `nct_ids` array
- [ ] Endpoint supports optional `state_filter`
- [ ] Returns refreshed trial data
- [ ] Handles API errors gracefully
- [ ] Returns empty data on failure (doesn't crash)

---

## **Error Handling**

- Invalid `nct_ids` → 400 Bad Request
- API failure → Return empty `trial_data` (graceful degradation)
- Timeout → Return partial results if any succeeded
- Log errors but don't fail frontend

---

**ESTIMATED TIME:** 20 minutes

