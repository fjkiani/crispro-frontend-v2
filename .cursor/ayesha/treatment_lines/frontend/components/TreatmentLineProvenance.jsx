/**
 * Treatment Line Provenance Display Component
 * 
 * Displays treatment line calculations and confidence adjustments
 * Shows line appropriateness, cross-resistance risk, sequencing fitness
 */

import React from 'react';
import PropTypes from 'prop-types';

const TreatmentLineProvenance = ({ provenance }) => {
    if (!provenance || provenance.error) {
        return null;
    }

    const {
        current_line,
        prior_therapies = [],
        line_appropriateness = 0,
        cross_resistance_risk = 0,
        sequencing_fitness = 0,
        nccn_category,
        confidence_penalty = 0,
        rationale
    } = provenance;

    // Color coding for scores
    const getScoreColor = (score, inverted = false) => {
        if (inverted) {
            // For cross-resistance risk (lower is better)
            if (score >= 0.5) return '#d32f2f'; // High risk - red
            if (score >= 0.3) return '#f57c00'; // Moderate risk - orange
            return '#388e3c'; // Low risk - green
        } else {
            // For line appropriateness and sequencing fitness (higher is better)
            if (score >= 0.8) return '#388e3c'; // Excellent - green
            if (score >= 0.6) return '#f57c00'; // Fair - orange
            return '#d32f2f'; // Poor - red
        }
    };

    const getScoreLabel = (score, inverted = false) => {
        if (inverted) {
            if (score >= 0.5) return 'High';
            if (score >= 0.3) return 'Moderate';
            return 'Low';
        } else {
            if (score >= 0.8) return 'Excellent';
            if (score >= 0.6) return 'Fair';
            return 'Poor';
        }
    };

    const getNccnBadgeColor = (category) => {
        switch(category) {
            case '1': return '#4caf50'; // Green - strongest evidence
            case '2A': return '#8bc34a'; // Light green
            case '2B': return '#ffc107'; // Amber
            case '3': return '#ff9800'; // Orange
            default: return '#9e9e9e'; // Grey
        }
    };

    return (
        <div className="treatment-line-provenance" style={{
            border: '1px solid #e0e0e0',
            borderRadius: '8px',
            padding: '16px',
            backgroundColor: '#fafafa',
            marginTop: '16px'
        }}>
            <div style={{ 
                display: 'flex', 
                justifyContent: 'space-between', 
                alignItems: 'center',
                marginBottom: '16px'
            }}>
                <h4 style={{ 
                    margin: 0, 
                    fontSize: '15px', 
                    fontWeight: 600,
                    color: '#333'
                }}>
                    Treatment Line Analysis
                </h4>
                <div style={{ display: 'flex', gap: '8px', alignItems: 'center' }}>
                    <span style={{
                        padding: '4px 8px',
                        backgroundColor: '#e3f2fd',
                        borderRadius: '12px',
                        fontSize: '12px',
                        fontWeight: 500,
                        color: '#1976d2'
                    }}>
                        Line {current_line}
                    </span>
                    {nccn_category && nccn_category !== 'unknown' && (
                        <span style={{
                            padding: '4px 8px',
                            backgroundColor: getNccnBadgeColor(nccn_category),
                            borderRadius: '12px',
                            fontSize: '12px',
                            fontWeight: 500,
                            color: 'white'
                        }}>
                            NCCN {nccn_category}
                        </span>
                    )}
                </div>
            </div>

            {/* Prior Therapies */}
            {prior_therapies.length > 0 && (
                <div style={{ marginBottom: '16px' }}>
                    <div style={{ fontSize: '12px', color: '#666', marginBottom: '6px' }}>
                        Prior Therapies:
                    </div>
                    <div style={{ display: 'flex', flexWrap: 'wrap', gap: '6px' }}>
                        {prior_therapies.map((therapy, index) => (
                            <span
                                key={index}
                                style={{
                                    padding: '4px 10px',
                                    backgroundColor: '#fff',
                                    border: '1px solid #ddd',
                                    borderRadius: '12px',
                                    fontSize: '12px',
                                    color: '#555'
                                }}
                            >
                                {therapy}
                            </span>
                        ))}
                    </div>
                </div>
            )}

            {/* Score Grid */}
            <div style={{ 
                display: 'grid', 
                gridTemplateColumns: 'repeat(3, 1fr)', 
                gap: '12px',
                marginBottom: '16px'
            }}>
                {/* Line Appropriateness */}
                <div style={{
                    padding: '12px',
                    backgroundColor: 'white',
                    borderRadius: '6px',
                    border: '1px solid #e0e0e0'
                }}>
                    <div style={{ fontSize: '11px', color: '#666', marginBottom: '6px' }}>
                        Line Fit
                    </div>
                    <div style={{ 
                        fontSize: '24px', 
                        fontWeight: 600,
                        color: getScoreColor(line_appropriateness)
                    }}>
                        {(line_appropriateness * 100).toFixed(0)}%
                    </div>
                    <div style={{ 
                        fontSize: '11px', 
                        color: getScoreColor(line_appropriateness),
                        fontWeight: 500,
                        marginTop: '4px'
                    }}>
                        {getScoreLabel(line_appropriateness)}
                    </div>
                </div>

                {/* Cross-Resistance Risk */}
                <div style={{
                    padding: '12px',
                    backgroundColor: 'white',
                    borderRadius: '6px',
                    border: '1px solid #e0e0e0'
                }}>
                    <div style={{ fontSize: '11px', color: '#666', marginBottom: '6px' }}>
                        Resistance Risk
                    </div>
                    <div style={{ 
                        fontSize: '24px', 
                        fontWeight: 600,
                        color: getScoreColor(cross_resistance_risk, true)
                    }}>
                        {(cross_resistance_risk * 100).toFixed(0)}%
                    </div>
                    <div style={{ 
                        fontSize: '11px', 
                        color: getScoreColor(cross_resistance_risk, true),
                        fontWeight: 500,
                        marginTop: '4px'
                    }}>
                        {getScoreLabel(cross_resistance_risk, true)} risk
                    </div>
                </div>

                {/* Sequencing Fitness */}
                <div style={{
                    padding: '12px',
                    backgroundColor: 'white',
                    borderRadius: '6px',
                    border: '1px solid #e0e0e0'
                }}>
                    <div style={{ fontSize: '11px', color: '#666', marginBottom: '6px' }}>
                        Sequencing Score
                    </div>
                    <div style={{ 
                        fontSize: '24px', 
                        fontWeight: 600,
                        color: getScoreColor(sequencing_fitness)
                    }}>
                        {(sequencing_fitness * 100).toFixed(0)}%
                    </div>
                    <div style={{ 
                        fontSize: '11px', 
                        color: getScoreColor(sequencing_fitness),
                        fontWeight: 500,
                        marginTop: '4px'
                    }}>
                        {getScoreLabel(sequencing_fitness)}
                    </div>
                </div>
            </div>

            {/* Confidence Penalty */}
            {confidence_penalty !== 0 && (
                <div style={{
                    padding: '12px',
                    backgroundColor: confidence_penalty > 0 ? '#fff3e0' : '#e8f5e9',
                    borderRadius: '6px',
                    border: '1px solid ' + (confidence_penalty > 0 ? '#ffb74d' : '#81c784'),
                    marginBottom: '12px'
                }}>
                    <div style={{ 
                        display: 'flex', 
                        alignItems: 'center',
                        gap: '8px'
                    }}>
                        <span style={{ fontSize: '16px' }}>
                            {confidence_penalty > 0 ? '⚠️' : '✓'}
                        </span>
                        <div style={{ flex: 1 }}>
                            <div style={{ 
                                fontSize: '13px', 
                                fontWeight: 500,
                                color: confidence_penalty > 0 ? '#e65100' : '#2e7d32'
                            }}>
                                Confidence {confidence_penalty > 0 ? 'Penalty' : 'Adjustment'}
                            </div>
                            <div style={{ fontSize: '12px', color: '#666', marginTop: '2px' }}>
                                {confidence_penalty > 0 ? '-' : '+'}{(Math.abs(confidence_penalty) * 100).toFixed(1)}%
                            </div>
                        </div>
                    </div>
                </div>
            )}

            {/* Rationale */}
            {rationale && (
                <div style={{
                    padding: '12px',
                    backgroundColor: 'white',
                    borderRadius: '6px',
                    border: '1px solid #e0e0e0'
                }}>
                    <div style={{ fontSize: '11px', color: '#666', marginBottom: '6px' }}>
                        Rationale:
                    </div>
                    <div style={{ fontSize: '13px', color: '#333', lineHeight: '1.5' }}>
                        {rationale}
                    </div>
                </div>
            )}

            {/* Tooltip explanations */}
            <div style={{
                marginTop: '12px',
                paddingTop: '12px',
                borderTop: '1px solid #e0e0e0',
                fontSize: '11px',
                color: '#666'
            }}>
                <strong>Line Fit:</strong> How appropriate is this drug for the current treatment line<br/>
                <strong>Resistance Risk:</strong> Risk of cross-resistance with prior therapies<br/>
                <strong>Sequencing Score:</strong> Overall sequencing fitness (combines fit and resistance)
            </div>
        </div>
    );
};

TreatmentLineProvenance.propTypes = {
    provenance: PropTypes.shape({
        current_line: PropTypes.number,
        prior_therapies: PropTypes.arrayOf(PropTypes.string),
        line_appropriateness: PropTypes.number,
        cross_resistance_risk: PropTypes.number,
        sequencing_fitness: PropTypes.number,
        nccn_category: PropTypes.string,
        confidence_penalty: PropTypes.number,
        rationale: PropTypes.string,
        error: PropTypes.string
    })
};

export default TreatmentLineProvenance;


