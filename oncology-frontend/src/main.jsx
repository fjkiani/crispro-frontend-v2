import React from "react";
import ReactDOM from "react-dom/client";
import { BrowserRouter as Router } from "react-router-dom";

import { StateContextProvider } from "./context";
import App from "./App";
import "./index.css";

// Global safety wrapper for Response.headers to prevent errors
// Some libraries may try to access headers on undefined or invalid responses
// This must run before any modules are imported to catch initialization errors

// 1. Protect Response.prototype.headers access
if (typeof Response !== 'undefined' && Response.prototype) {
  const originalHeadersDescriptor = Object.getOwnPropertyDescriptor(Response.prototype, 'headers');
  
  // Ensure headers getter always returns a valid Headers object
  Object.defineProperty(Response.prototype, 'headers', {
    get: function() {
      // If this is not a valid Response instance, return empty Headers
      if (!this || !(this instanceof Response)) {
        console.warn('Accessing headers on invalid Response object');
        return new Headers();
      }
      // Return original headers if available
      if (originalHeadersDescriptor && originalHeadersDescriptor.get) {
        try {
          return originalHeadersDescriptor.get.call(this);
        } catch (e) {
          console.warn('Error accessing headers:', e);
          return new Headers();
        }
      }
      return new Headers();
    },
    configurable: true,
    enumerable: true
  });
}

// 2. Wrap fetch to ensure response objects are always valid
if (typeof window !== 'undefined' && window.fetch) {
  const originalFetch = window.fetch;
  window.fetch = async function(...args) {
    try {
      const response = await originalFetch.apply(this, args);
      // Ensure response is valid and has headers
      if (!response) {
        throw new Error('Fetch returned undefined response');
      }
      // Headers should always exist, but ensure it does
      if (!response.headers) {
        console.warn('Response missing headers, creating empty Headers');
        Object.defineProperty(response, 'headers', {
          value: new Headers(),
          writable: false,
          configurable: true
        });
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
