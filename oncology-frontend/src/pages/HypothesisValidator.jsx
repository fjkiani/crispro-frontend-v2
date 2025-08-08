import React, { useState } from 'react';
import { Box, TextField, Button, Typography, CircularProgress, Paper, Link } from '@mui/material';

const HypothesisValidator = () => {
    const [query, setQuery] = useState('What is the mechanism of action for Neovastat (AE-941)?');
    const [isLoading, setIsLoading] = useState(false);
    const [results, setResults] = useState(null);
    const [error, setError] = useState(null);

    const handleSearch = async () => {
        setIsLoading(true);
        setError(null);
        setResults(null);

        try {
            const response = await fetch('http://localhost:8000/api/intelligence/search', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({ query }),
            });

            if (!response.ok) {
                throw new Error(`API request failed with status ${response.status}`);
            }

            const data = await response.json();
            setResults(data);
        } catch (err) {
            setError(err.message);
        } finally {
            setIsLoading(false);
        }
    };

    return (
        <Box sx={{ p: 4 }}>
            <Typography variant="h4" gutterBottom>
                🔬 Hypothesis Validator
            </Typography>
            <Typography variant="subtitle1" color="text.secondary" gutterBottom>
                OPERATION: DEEP DIVE - Manual Reconnaissance (Zeta-Stream Ingress Protocol)
            </Typography>

            <Paper sx={{ p: 3, mt: 2 }}>
                <Typography variant="h6">New Intelligence Run</Typography>
                <Box sx={{ display: 'flex', alignItems: 'center', mt: 2, gap: 2 }}>
                    <TextField
                        fullWidth
                        label="Research Query"
                        variant="outlined"
                        value={query}
                        onChange={(e) => setQuery(e.target.value)}
                        disabled={isLoading}
                    />
                    <Button
                        variant="contained"
                        onClick={handleSearch}
                        disabled={isLoading}
                        size="large"
                    >
                        {isLoading ? <CircularProgress size={24} /> : 'Execute'}
                    </Button>
                </Box>
            </Paper>

            {error && (
                <Paper sx={{ p: 3, mt: 4, bgcolor: 'error.main', color: 'error.contrastText' }}>
                    <Typography variant="h6">Error</Typography>
                    <Typography>{error}</Typography>
                </Paper>
            )}

            {results && (
                <Box sx={{ mt: 4 }}>
                    {results.answer && (
                        <Paper sx={{ p: 3, mb: 4, bgcolor: 'info.main', color: 'info.contrastText' }}>
                            <Typography variant="h6" gutterBottom>Synthesized Answer</Typography>
                            <Typography>{results.answer}</Typography>
                        </Paper>
                    )}

                    <Typography variant="h5" gutterBottom>Primary Sources</Typography>
                    {results.results && results.results.map((item, index) => (
                        <Paper key={index} sx={{ p: 3, mb: 2 }}>
                            <Typography variant="h6">
                                <Link href={item.url} target="_blank" rel="noopener">
                                    {item.title}
                                </Link>
                            </Typography>
                            <Typography variant="caption" color="text.secondary" gutterBottom>
                                Source: {item.url} | Score: {item.score.toFixed(2)}
                            </Typography>
                            <Typography variant="body1" sx={{ mt: 1 }}>
                                {item.content}
                            </Typography>
                        </Paper>
                    ))}
                </Box>
            )}
        </Box>
    );
};

export default HypothesisValidator; 