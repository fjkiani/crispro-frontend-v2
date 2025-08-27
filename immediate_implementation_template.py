# Immediate Implementation Template for Data Source Improvements

import asyncio
import time
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
import json
import logging
from datetime import datetime

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class DataSourceStatus:
    name: str
    status: str  # 'healthy', 'degraded', 'down'
    last_check: datetime
    response_time: float
    error_count: int
    consecutive_failures: int

class CircuitBreaker:
    def __init__(self, failure_threshold: int = 5, recovery_timeout: int = 60):
        self.failure_threshold = failure_threshold
        self.recovery_timeout = recovery_timeout
        self.failure_count = 0
        self.last_failure_time = None
        self.state = 'CLOSED'  # CLOSED, OPEN, HALF_OPEN
    
    def is_available(self) -> bool:
        if self.state == 'CLOSED':
            return True
        elif self.state == 'OPEN':
            if self.last_failure_time and \
               (time.time() - self.last_failure_time) > self.recovery_timeout:
                self.state = 'HALF_OPEN'
                return True
            return False
        elif self.state == 'HALF_OPEN':
            return True
        return False
    
    def record_success(self):
        self.failure_count = 0
        self.state = 'CLOSED'
    
    def record_failure(self):
        self.failure_count += 1
        self.last_failure_time = time.time()
        
        if self.failure_count >= self.failure_threshold:
            self.state = 'OPEN'

class AdaptiveRateLimiter:
    def __init__(self):
        self.requests_per_minute = {
            'ncbi': 10,
            'ensembl': 15,
            'clinvar': 5,
            'pubmed': 10,
            'cosmic': 20
        }
        self.current_usage = {}
        self.backoff_multipliers = {}
    
    async def wait_if_needed(self, source: str):
        now = time.time()
        window_start = now - 60
        
        if source not in self.current_usage:
            self.current_usage[source] = []
        
        # Clean old requests
        self.current_usage[source] = [
            t for t in self.current_usage[source] if t > window_start
        ]
        
        current_usage = len(self.current_usage[source])
        limit = self.requests_per_minute[source]
        
        if current_usage >= limit:
            # Calculate wait time with exponential backoff
            backoff = self.backoff_multipliers.get(source, 1)
            wait_time = 60 + (backoff * 10)  # Base 60s + backoff
            
            logger.warning(f"Rate limit hit for {source}, waiting {wait_time}s")
            await asyncio.sleep(wait_time)
            
            # Increase backoff for next time
            self.backoff_multipliers[source] = min(backoff * 2, 10)
        else:
            # Reset backoff on successful request window
            self.backoff_multipliers[source] = 1
        
        self.current_usage[source].append(now)

class UnifiedDataClient:
    """
    Immediate implementation template for improved data source management.
    Copy this into your existing codebase and extend it.
    """
    
    def __init__(self):
        self.circuit_breakers = {
            'ncbi': CircuitBreaker(),
            'clinvar': CircuitBreaker(),
            'ensembl': CircuitBreaker(),
            'pubmed': CircuitBreaker(),
            'cosmic': CircuitBreaker()
        }
        
        self.rate_limiter = AdaptiveRateLimiter()
        self.source_status = {}
        self.cache = {}
    
    async def get_gene_data(self, gene_symbol: str) -> Dict[str, Any]:
        """
        Example method showing improved data retrieval pattern.
        Replace your existing NCBI client calls with this pattern.
        """
        results = {}
        
        # Try NCBI first (your current working source)
        try:
            if self.circuit_breakers['ncbi'].is_available():
                await self.rate_limiter.wait_if_needed('ncbi')
                
                # Replace with your actual NCBI client call
                ncbi_data = await self._fetch_from_ncbi(gene_symbol)
                
                if ncbi_data:
                    results['ncbi'] = ncbi_data
                    self.circuit_breakers['ncbi'].record_success()
                else:
                    self.circuit_breakers['ncbi'].record_failure()
                    
        except Exception as e:
            logger.error(f"NCBI failed for {gene_symbol}: {e}")
            self.circuit_breakers['ncbi'].record_failure()
        
        # Try Ensembl as fallback
        if not results and self.circuit_breakers['ensembl'].is_available():
            try:
                await self.rate_limiter.wait_if_needed('ensembl')
                
                # TODO: Implement Ensembl client
                ensembl_data = await self._fetch_from_ensembl(gene_symbol)
                
                if ensembl_data:
                    results['ensembl'] = ensembl_data
                    self.circuit_breakers['ensembl'].record_success()
                    
            except Exception as e:
                logger.error(f"Ensembl failed for {gene_symbol}: {e}")
                self.circuit_breakers['ensembl'].record_failure()
        
        # Try COSMIC for cancer-specific data
        if self.circuit_breakers['cosmic'].is_available():
            try:
                await self.rate_limiter.wait_if_needed('cosmic')
                
                # Use your existing COSMIC client
                cosmic_data = await self._fetch_from_cosmic(gene_symbol)
                
                if cosmic_data:
                    results['cosmic'] = cosmic_data
                    self.circuit_breakers['cosmic'].record_success()
                    
            except Exception as e:
                logger.error(f"COSMIC failed for {gene_symbol}: {e}")
                self.circuit_breakers['cosmic'].record_failure()
        
        # Merge and return results
        return self._merge_gene_data(results)
    
    async def _fetch_from_ncbi(self, gene_symbol: str) -> Optional[Dict]:
        """Replace with your existing NCBI client logic"""
        # This is where you put your current NCBI client code
        # Return None if no data found
        return None
    
    async def _fetch_from_ensembl(self, gene_symbol: str) -> Optional[Dict]:
        """Implement Ensembl client - placeholder for future"""
        # TODO: Implement Ensembl REST API client
        # https://rest.ensembl.org/documentation/info/symbol_lookup
        return None
    
    async def _fetch_from_cosmic(self, gene_symbol: str) -> Optional[Dict]:
        """Use your existing COSMIC client"""
        # This is where you put your current COSMIC client code
        return None
    
    def _merge_gene_data(self, sources: Dict[str, Dict]) -> Dict[str, Any]:
        """Merge data from multiple sources into unified format"""
        merged = {
            'gene_symbol': None,
            'sources': list(sources.keys()),
            'data': {},
            'last_updated': datetime.now().isoformat()
        }
        
        # Extract gene symbol from any source
        for source_data in sources.values():
            if 'gene_symbol' in source_data:
                merged['gene_symbol'] = source_data['gene_symbol']
                break
        
        # Add source-specific data
        merged['data'] = sources
        
        return merged
    
    async def get_health_status(self) -> Dict[str, DataSourceStatus]:
        """Get health status of all data sources"""
        status = {}
        for source_name, breaker in self.circuit_breakers.items():
            status[source_name] = DataSourceStatus(
                name=source_name,
                status='healthy' if breaker.state == 'CLOSED' else 'degraded',
                last_check=datetime.now(),
                response_time=0.0,  # TODO: Track actual response times
                error_count=breaker.failure_count,
                consecutive_failures=breaker.failure_count
            )
        
        return status

# Example usage in your existing code:
async def example_usage():
    """How to integrate this into your existing workflow"""
    client = UnifiedDataClient()
    
    # Instead of calling NCBI directly:
    # old_way = ncbi_client.find_protein_domains(gene_symbol)
    
    # New way with resilience and fallbacks:
    gene_data = await client.get_gene_data("BRCA1")
    
    # The result includes data from all available sources
    if gene_data['sources']:
        print(f"Found data from: {gene_data['sources']}")
        # Use the merged data in your existing workflow
        domains = gene_data['data'].get('ncbi', {}).get('domains', [])
        
        # Your existing logic continues here...
    
    # Monitor health of data sources
    health = await client.get_health_status()
    for source, status in health.items():
        if status.status == 'degraded':
            logger.warning(f"Data source {source} is degraded")

if __name__ == "__main__":
    # Quick test
    asyncio.run(example_usage())
