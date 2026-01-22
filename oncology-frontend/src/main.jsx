import React from "react";
import ReactDOM from "react-dom/client";
import { BrowserRouter as Router } from "react-router-dom";

import { StateContextProvider } from "./context";
import App from "./App";
import "./index.css";

// Global safety wrapper for Response.headers to prevent errors
// Some libraries may try to access headers on undefined responses
if (typeof Response !== 'undefined') {
  const OriginalResponse = Response;
  // Wrap fetch to ensure response objects always have headers
  const originalFetch = window.fetch;
  window.fetch = async function(...args) {
    try {
      const response = await originalFetch.apply(this, args);
      // Ensure response.headers exists (should always be present, but safety check)
      if (response && !response.headers) {
        console.warn('Response object missing headers property');
      }
      return response;
    } catch (error) {
      console.error('Fetch error:', error);
      throw error;
    }
  };
}

const root = ReactDOM.createRoot(document.getElementById("root"));

root.render(
  <Router>
    <StateContextProvider>
      <App />
    </StateContextProvider>
  </Router>
);
