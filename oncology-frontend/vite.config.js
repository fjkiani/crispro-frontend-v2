import path from "path";
import react from "@vitejs/plugin-react";
import { defineConfig } from "vite";
import { Buffer } from 'buffer';

export default defineConfig({
  plugins: [react()],
  logLevel: 'warn', // Reduce build output noise (suppresses 'use client' warnings)
  resolve: {
    alias: {
      "@": path.resolve(__dirname, "./src"),
      'buffer': 'buffer',
      'path': 'path-browserify',
      '~': path.resolve(__dirname, "./src"),
      '#minpath': path.resolve(__dirname, "./node_modules/vfile/lib/minpath.js"),
      '#minproc': path.resolve(__dirname, "./node_modules/vfile/lib/minproc.js"),
      '#minurl': path.resolve(__dirname, "./node_modules/vfile/lib/minurl.js"),
      'node:url': 'url',
    },
  },
  define: {
    'global.Buffer': Buffer,
    global: "globalThis",
    "process.env": {},
  },
  base: "/",
  server: {
    // No proxy needed for standalone frontend deployment
  },
  build: {
    // Increase chunk size warning limit
    chunkSizeWarningLimit: 1000,
    rollupOptions: {
      external: [],
      output: {
        // Manual chunk splitting to reduce memory usage
        // Conservative approach: only split very large, independent libraries
        manualChunks: (id) => {
          if (id.includes('node_modules')) {
            // Keep MUI packages together (they have tight dependencies)
            if (id.includes('@mui/material') || id.includes('@mui/icons-material') || 
                id.includes('@mui/system') || id.includes('@mui/lab') || 
                id.includes('@mui/x-charts') || id.includes('@emotion')) {
              return 'mui';
            }
            // React core libraries together
            if (id.includes('react') || id.includes('react-dom') || id.includes('react-router')) {
              return 'react-vendor';
            }
            // Supabase (independent)
            if (id.includes('@supabase')) {
              return 'supabase';
            }
            // Recharts (independent)
            if (id.includes('recharts')) {
              return 'recharts';
            }
            // AI libraries (independent)
            if (id.includes('@google/generative-ai') || id.includes('openai')) {
              return 'ai-vendor';
            }
            // Web3 libraries (independent)
            if (id.includes('ethers') || id.includes('@thirdweb')) {
              return 'web3-vendor';
            }
            // Everything else goes to vendor (less aggressive splitting)
            return 'vendor';
          }
        },
      },
    },
    // Reduce memory usage during build
    minify: 'esbuild', // Faster than terser
    target: 'es2015', // Broader compatibility, less transformation
  },
});
