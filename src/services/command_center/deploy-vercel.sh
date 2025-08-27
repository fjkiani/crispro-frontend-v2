#!/bin/bash

# Command Center Backend Vercel Deployment Script
echo "🚀 Deploying Command Center Backend to Vercel..."

# Check if we're in the right directory
if [ ! -f "vercel_main.py" ]; then
    echo "❌ vercel_main.py not found. Please run this script from the command_center directory."
    exit 1
fi

# Check if vercel CLI is installed
if ! command -v vercel &> /dev/null; then
    echo "❌ Vercel CLI not found. Installing..."
    npm install -g vercel
fi

# Check if requirements-vercel.txt exists
if [ ! -f "requirements-vercel.txt" ]; then
    echo "❌ requirements-vercel.txt not found. Please ensure it exists."
    exit 1
fi

# Deploy to Vercel
echo "📦 Deploying to Vercel..."
vercel --prod

echo "✅ Command Center Backend deployment to Vercel complete!"
echo "🔗 Your backend will be available at the Vercel URL shown above"
echo "📝 Update your frontend .env file with the backend Vercel URL"
echo "💡 Note: Vercel has a 10-second timeout limit for serverless functions"
echo "⚠️  Warning: This service may hit timeout limits for long-running AI operations" 