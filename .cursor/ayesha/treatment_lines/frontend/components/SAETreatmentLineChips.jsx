/**
 * SAE Treatment Line Chips Component
 * 
 * Displays the 3 treatment line SAE features:
 * 1. Line Appropriateness (line_appropriateness)
 * 2. Cross-Resistance Risk (cross_resistance_risk)
 * 3. Sequencing Fitness (sequencing_fitness)
 */

import React from 'react';
import PropTypes from 'prop-types';

const SAETreatmentLineChips = ({ saeFeatures }) => {
    if (!saeFeatures || saeFeatures.length === 0) {
        return null;
    }

    // Filter to treatment line features
    const treatmentLineFeatures = saeFeatures.filter(feature => 
        ['line_appropriateness', 'cross_resistance_risk', 'sequencing_fitness'].includes(feature.id)
    );

    if (treatmentLineFeatures.length === 0) {
        return null;
    }

    const getChipColor = (feature) => {
        const { impact, activation } = feature;
        
        // For cross_resistance_risk (inverted - high is bad)
        if (feature.id === 'cross_resistance_risk') {
            if (activation >= 0.5) return { bg: '#ffebee', border: '#ef5350', text: '#c62828' }; // High risk - red
            if (activation >= 0.3) return { bg: '#fff3e0', border: '#ff9800', text: '#e65100' }; // Moderate - orange
            return { bg: '#e8f5e9', border: '#66bb6a', text: '#2e7d32' }; // Low risk - green
        }
        
        // For line_appropriateness and sequencing_fitness (higher is better)
        if (impact === 'positive') {
            return { bg: '#e8f5e9', border: '#66bb6a', text: '#2e7d32' }; // Green
        } else if (impact === 'negative') {
            return { bg: '#ffebee', border: '#ef5350', text: '#c62828' }; // Red
        } else {
            return { bg: '#fff3e0', border: '#ff9800', text: '#e65100' }; // Orange (neutral/warning)
        }
    };

    const getIcon = (feature) => {
        switch(feature.id) {
            case 'line_appropriateness':
                return feature.impact === 'positive' ? '✓' : '⚠';
            case 'cross_resistance_risk':
                return feature.activation >= 0.5 ? '⚠️' : feature.activation >= 0.3 ? '⚡' : '✓';
            case 'sequencing_fitness':
                return feature.impact === 'positive' ? '⭐' : '⚠';
            default:
                return '•';
        }
    };

    const getLabel = (feature) => {
        switch(feature.id) {
            case 'line_appropriateness':
                return 'Line Fit';
            case 'cross_resistance_risk':
                return 'Resistance Risk';
            case 'sequencing_fitness':
                return 'Sequencing Score';
            default:
                return feature.name || feature.id;
        }
    };

    return (
        <div className="sae-treatment-line-chips" style={{
            display: 'flex',
            flexWrap: 'wrap',
            gap: '8px',
            marginTop: '12px'
        }}>
            <div style={{
                width: '100%',
                fontSize: '12px',
                fontWeight: 500,
                color: '#666',
                marginBottom: '4px'
            }}>
                Treatment Line Features:
            </div>
            
            {treatmentLineFeatures.map((feature, index) => {
                const colors = getChipColor(feature);
                const icon = getIcon(feature);
                const label = getLabel(feature);
                const value = (feature.activation * 100).toFixed(0);
                
                return (
                    <div
                        key={index}
                        style={{
                            display: 'inline-flex',
                            alignItems: 'center',
                            padding: '8px 14px',
                            backgroundColor: colors.bg,
                            border: `1px solid ${colors.border}`,
                            borderRadius: '20px',
                            fontSize: '13px',
                            fontWeight: 500,
                            color: colors.text,
                            cursor: 'pointer',
                            transition: 'all 0.2s ease'
                        }}
                        title={feature.explanation || ''}
                    >
                        <span style={{ marginRight: '6px', fontSize: '14px' }}>
                            {icon}
                        </span>
                        <span>
                            {label}: {value}%
                        </span>
                    </div>
                );
            })}

            {/* Explanation tooltip */}
            <div style={{
                width: '100%',
                fontSize: '11px',
                color: '#999',
                marginTop: '4px',
                lineHeight: '1.4'
            }}>
                Hover over chips to see detailed explanations
            </div>
        </div>
    );
};

SAETreatmentLineChips.propTypes = {
    saeFeatures: PropTypes.arrayOf(PropTypes.shape({
        id: PropTypes.string.isRequired,
        name: PropTypes.string,
        activation: PropTypes.number.isRequired,
        impact: PropTypes.string.isRequired,
        explanation: PropTypes.string
    }))
};

export default SAETreatmentLineChips;


