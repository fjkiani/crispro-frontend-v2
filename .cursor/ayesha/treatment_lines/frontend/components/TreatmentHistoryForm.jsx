/**
 * Treatment History Form Component
 * 
 * Collects patient treatment line context for efficacy prediction.
 * Sends treatment_history to backend efficacy endpoint.
 */

import React, { useState } from 'react';
import PropTypes from 'prop-types';

const TreatmentHistoryForm = ({ 
    onSubmit, 
    defaultLine = 1, 
    defaultPriorTherapies = [],
    disease = "unknown"
}) => {
    const [currentLine, setCurrentLine] = useState(defaultLine);
    const [priorTherapies, setPriorTherapies] = useState(defaultPriorTherapies);
    const [newTherapy, setNewTherapy] = useState('');
    const [showForm, setShowForm] = useState(false);

    // Common therapies by disease for quick add
    const therapySuggestions = {
        ovarian_cancer: [
            'carboplatin',
            'paclitaxel',
            'carboplatin+paclitaxel',
            'bevacizumab',
            'olaparib',
            'niraparib',
            'rucaparib'
        ],
        breast_her2_positive: [
            'trastuzumab',
            'pertuzumab',
            'paclitaxel',
            'docetaxel',
            'trastuzumab deruxtecan',
            'tucatinib',
            'capecitabine'
        ],
        breast_triple_negative: [
            'carboplatin',
            'paclitaxel',
            'pembrolizumab',
            'sacituzumab govitecan',
            'capecitabine'
        ]
    };

    const suggestions = therapySuggestions[disease] || [];

    const handleAddTherapy = (therapy) => {
        if (therapy && !priorTherapies.includes(therapy)) {
            setPriorTherapies([...priorTherapies, therapy]);
            setNewTherapy('');
        }
    };

    const handleRemoveTherapy = (index) => {
        setPriorTherapies(priorTherapies.filter((_, i) => i !== index));
    };

    const handleSubmit = (e) => {
        e.preventDefault();
        
        const treatmentHistory = {
            current_line: currentLine,
            prior_therapies: priorTherapies
        };

        onSubmit(treatmentHistory);
    };

    const handleReset = () => {
        setCurrentLine(1);
        setPriorTherapies([]);
        setNewTherapy('');
    };

    if (!showForm) {
        return (
            <div className="treatment-history-toggle">
                <button
                    onClick={() => setShowForm(true)}
                    className="btn-secondary"
                    style={{
                        padding: '8px 16px',
                        backgroundColor: '#f0f0f0',
                        border: '1px solid #ccc',
                        borderRadius: '4px',
                        cursor: 'pointer',
                        fontSize: '14px'
                    }}
                >
                    + Add Treatment History (Line {currentLine})
                </button>
                {priorTherapies.length > 0 && (
                    <span style={{ marginLeft: '12px', fontSize: '12px', color: '#666' }}>
                        {priorTherapies.length} prior therap{priorTherapies.length === 1 ? 'y' : 'ies'}
                    </span>
                )}
            </div>
        );
    }

    return (
        <div className="treatment-history-form" style={{
            border: '1px solid #ddd',
            borderRadius: '8px',
            padding: '20px',
            backgroundColor: '#fafafa',
            marginBottom: '20px'
        }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '16px' }}>
                <h3 style={{ margin: 0, fontSize: '16px', fontWeight: 600 }}>
                    Treatment History
                </h3>
                <button
                    onClick={() => setShowForm(false)}
                    style={{
                        background: 'none',
                        border: 'none',
                        cursor: 'pointer',
                        fontSize: '18px',
                        color: '#666'
                    }}
                >
                    ✕
                </button>
            </div>

            <form onSubmit={handleSubmit}>
                {/* Current Treatment Line */}
                <div style={{ marginBottom: '20px' }}>
                    <label style={{ display: 'block', marginBottom: '8px', fontWeight: 500, fontSize: '14px' }}>
                        Current Treatment Line *
                    </label>
                    <select
                        value={currentLine}
                        onChange={(e) => setCurrentLine(parseInt(e.target.value))}
                        required
                        style={{
                            width: '100%',
                            padding: '8px 12px',
                            border: '1px solid #ccc',
                            borderRadius: '4px',
                            fontSize: '14px'
                        }}
                    >
                        <option value={1}>1st Line (First-line therapy)</option>
                        <option value={2}>2nd Line (Post-progression)</option>
                        <option value={3}>3rd Line (Second progression)</option>
                        <option value={4}>4th Line (Third progression)</option>
                        <option value={5}>5th Line (Fourth progression)</option>
                        <option value={6}>6th Line</option>
                        <option value={7}>7th Line</option>
                        <option value={8}>8th Line</option>
                        <option value={9}>9th Line</option>
                        <option value={10}>10th Line</option>
                    </select>
                    <div style={{ fontSize: '12px', color: '#666', marginTop: '4px' }}>
                        Select the line of therapy you are planning
                    </div>
                </div>

                {/* Prior Therapies */}
                <div style={{ marginBottom: '20px' }}>
                    <label style={{ display: 'block', marginBottom: '8px', fontWeight: 500, fontSize: '14px' }}>
                        Prior Therapies {currentLine > 1 && '*'}
                    </label>
                    
                    {/* Current prior therapies list */}
                    {priorTherapies.length > 0 && (
                        <div style={{ marginBottom: '12px' }}>
                            {priorTherapies.map((therapy, index) => (
                                <div
                                    key={index}
                                    style={{
                                        display: 'inline-flex',
                                        alignItems: 'center',
                                        padding: '6px 12px',
                                        margin: '4px',
                                        backgroundColor: '#e3f2fd',
                                        borderRadius: '16px',
                                        fontSize: '13px'
                                    }}
                                >
                                    <span>{therapy}</span>
                                    <button
                                        type="button"
                                        onClick={() => handleRemoveTherapy(index)}
                                        style={{
                                            marginLeft: '8px',
                                            background: 'none',
                                            border: 'none',
                                            cursor: 'pointer',
                                            fontSize: '16px',
                                            color: '#666'
                                        }}
                                    >
                                        ✕
                                    </button>
                                </div>
                            ))}
                        </div>
                    )}

                    {/* Add new therapy input */}
                    <div style={{ display: 'flex', gap: '8px', marginBottom: '12px' }}>
                        <input
                            type="text"
                            value={newTherapy}
                            onChange={(e) => setNewTherapy(e.target.value)}
                            onKeyPress={(e) => {
                                if (e.key === 'Enter') {
                                    e.preventDefault();
                                    handleAddTherapy(newTherapy);
                                }
                            }}
                            placeholder="Enter drug or regimen name"
                            style={{
                                flex: 1,
                                padding: '8px 12px',
                                border: '1px solid #ccc',
                                borderRadius: '4px',
                                fontSize: '14px'
                            }}
                        />
                        <button
                            type="button"
                            onClick={() => handleAddTherapy(newTherapy)}
                            disabled={!newTherapy}
                            style={{
                                padding: '8px 16px',
                                backgroundColor: newTherapy ? '#2196f3' : '#ccc',
                                color: 'white',
                                border: 'none',
                                borderRadius: '4px',
                                cursor: newTherapy ? 'pointer' : 'not-allowed',
                                fontSize: '14px'
                            }}
                        >
                            Add
                        </button>
                    </div>

                    {/* Quick add suggestions */}
                    {suggestions.length > 0 && (
                        <div>
                            <div style={{ fontSize: '12px', color: '#666', marginBottom: '6px' }}>
                                Quick add:
                            </div>
                            <div style={{ display: 'flex', flexWrap: 'wrap', gap: '6px' }}>
                                {suggestions
                                    .filter(s => !priorTherapies.includes(s))
                                    .slice(0, 5)
                                    .map((suggestion, index) => (
                                        <button
                                            key={index}
                                            type="button"
                                            onClick={() => handleAddTherapy(suggestion)}
                                            style={{
                                                padding: '4px 10px',
                                                backgroundColor: 'white',
                                                border: '1px solid #2196f3',
                                                borderRadius: '12px',
                                                cursor: 'pointer',
                                                fontSize: '12px',
                                                color: '#2196f3'
                                            }}
                                        >
                                            + {suggestion}
                                        </button>
                                    ))
                                }
                            </div>
                        </div>
                    )}

                    <div style={{ fontSize: '12px', color: '#666', marginTop: '8px' }}>
                        {currentLine === 1 
                            ? 'First-line therapy - no prior therapies needed'
                            : `Add all prior therapies received before line ${currentLine}`
                        }
                    </div>
                </div>

                {/* Action buttons */}
                <div style={{ display: 'flex', gap: '12px', marginTop: '20px' }}>
                    <button
                        type="submit"
                        disabled={currentLine > 1 && priorTherapies.length === 0}
                        style={{
                            flex: 1,
                            padding: '10px',
                            backgroundColor: (currentLine > 1 && priorTherapies.length === 0) ? '#ccc' : '#4caf50',
                            color: 'white',
                            border: 'none',
                            borderRadius: '4px',
                            cursor: (currentLine > 1 && priorTherapies.length === 0) ? 'not-allowed' : 'pointer',
                            fontSize: '14px',
                            fontWeight: 500
                        }}
                    >
                        Apply Treatment History
                    </button>
                    <button
                        type="button"
                        onClick={handleReset}
                        style={{
                            padding: '10px 20px',
                            backgroundColor: 'white',
                            color: '#666',
                            border: '1px solid #ccc',
                            borderRadius: '4px',
                            cursor: 'pointer',
                            fontSize: '14px'
                        }}
                    >
                        Reset
                    </button>
                </div>
            </form>

            {/* Validation warning */}
            {currentLine > 1 && priorTherapies.length === 0 && (
                <div style={{
                    marginTop: '12px',
                    padding: '12px',
                    backgroundColor: '#fff3cd',
                    border: '1px solid #ffc107',
                    borderRadius: '4px',
                    fontSize: '13px',
                    color: '#856404'
                }}>
                    ⚠️ Please add at least one prior therapy for line {currentLine}
                </div>
            )}
        </div>
    );
};

TreatmentHistoryForm.propTypes = {
    onSubmit: PropTypes.func.isRequired,
    defaultLine: PropTypes.number,
    defaultPriorTherapies: PropTypes.arrayOf(PropTypes.string),
    disease: PropTypes.string
};

export default TreatmentHistoryForm;









