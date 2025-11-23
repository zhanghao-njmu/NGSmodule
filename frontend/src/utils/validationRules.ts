/**
 * Validation Rules Utility
 * Centralized form validation rules for consistent validation across the application
 * Compatible with Ant Design Form component
 */

import type { Rule } from 'antd/es/form'

/**
 * Validation rule library
 * Provides reusable validation rules for common form fields
 */
export const validationRules = {
  /**
   * Required field validation
   *
   * @param fieldName - Human-readable field name for error message
   * @returns Validation rule
   *
   * @example
   * <Form.Item name="username" rules={[validationRules.required('username')]}>
   */
  required: (fieldName: string): Rule => ({
    required: true,
    message: `Please input ${fieldName}!`,
    whitespace: true,
  }),

  /**
   * Email validation
   * Validates standard email format
   */
  email: {
    type: 'email' as const,
    message: 'Please input a valid email address!',
  },

  /**
   * URL validation
   * Validates standard URL format
   */
  url: {
    type: 'url' as const,
    message: 'Please input a valid URL!',
  },

  /**
   * Minimum length validation
   *
   * @param min - Minimum number of characters
   * @returns Validation rule
   *
   * @example
   * rules={[validationRules.minLength(8)]}
   */
  minLength: (min: number): Rule => ({
    min,
    message: `Must be at least ${min} characters!`,
  }),

  /**
   * Maximum length validation
   *
   * @param max - Maximum number of characters
   * @returns Validation rule
   */
  maxLength: (max: number): Rule => ({
    max,
    message: `Must be at most ${max} characters!`,
  }),

  /**
   * Exact length validation
   *
   * @param length - Required length
   * @returns Validation rule
   */
  exactLength: (length: number): Rule => ({
    len: length,
    message: `Must be exactly ${length} characters!`,
  }),

  /**
   * Length range validation
   *
   * @param min - Minimum length
   * @param max - Maximum length
   * @returns Validation rule
   */
  lengthRange: (min: number, max: number): Rule => ({
    min,
    max,
    message: `Must be between ${min} and ${max} characters!`,
  }),

  /**
   * Number minimum value validation
   *
   * @param min - Minimum value
   * @returns Validation rule
   */
  minValue: (min: number): Rule => ({
    type: 'number' as const,
    min,
    message: `Must be at least ${min}!`,
  }),

  /**
   * Number maximum value validation
   *
   * @param max - Maximum value
   * @returns Validation rule
   */
  maxValue: (max: number): Rule => ({
    type: 'number' as const,
    max,
    message: `Must be at most ${max}!`,
  }),

  /**
   * Number range validation
   *
   * @param min - Minimum value
   * @param max - Maximum value
   * @returns Validation rule
   */
  numberRange: (min: number, max: number): Rule => ({
    type: 'number' as const,
    min,
    max,
    message: `Must be between ${min} and ${max}!`,
  }),

  /**
   * Pattern matching validation
   *
   * @param pattern - Regular expression pattern
   * @param message - Custom error message
   * @returns Validation rule
   *
   * @example
   * validationRules.pattern(/^[A-Z]/, 'Must start with uppercase letter')
   */
  pattern: (pattern: RegExp, message: string): Rule => ({
    pattern,
    message,
  }),

  /**
   * Username validation rules
   * - Required
   * - 3-50 characters
   * - Alphanumeric, hyphens, and underscores only
   */
  username: [
    {
      required: true,
      message: 'Please input your username!',
      whitespace: true,
    },
    {
      min: 3,
      message: 'Username must be at least 3 characters!',
    },
    {
      max: 50,
      message: 'Username must be at most 50 characters!',
    },
    {
      pattern: /^[a-zA-Z0-9_-]+$/,
      message: 'Username can only contain letters, numbers, hyphens, and underscores!',
    },
  ],

  /**
   * Password validation rules
   * - Required
   * - Minimum 8 characters
   * - At least one uppercase letter
   * - At least one lowercase letter
   * - At least one number
   * - At least one special character (optional, commented out)
   */
  password: [
    {
      required: true,
      message: 'Please input your password!',
    },
    {
      min: 8,
      message: 'Password must be at least 8 characters!',
    },
    {
      pattern: /[A-Z]/,
      message: 'Password must contain at least one uppercase letter!',
    },
    {
      pattern: /[a-z]/,
      message: 'Password must contain at least one lowercase letter!',
    },
    {
      pattern: /[0-9]/,
      message: 'Password must contain at least one number!',
    },
    // Uncomment for special character requirement
    // {
    //   pattern: /[!@#$%^&*(),.?":{}|<>]/,
    //   message: 'Password must contain at least one special character!',
    // },
  ],

  /**
   * Simple password validation (less strict)
   * - Required
   * - Minimum 6 characters
   */
  passwordSimple: [
    {
      required: true,
      message: 'Please input your password!',
    },
    {
      min: 6,
      message: 'Password must be at least 6 characters!',
    },
  ],

  /**
   * Confirm password validation
   * Must match the 'password' field
   *
   * @param getFieldValue - Form's getFieldValue function
   * @param passwordField - Name of the password field to match (default: 'password')
   * @returns Validation rule
   *
   * @example
   * const form = Form.useForm()
   * <Form.Item
   *   name="confirmPassword"
   *   rules={[validationRules.confirmPassword(form.getFieldValue)]}
   * >
   */
  confirmPassword: (getFieldValue: (field: string) => any, passwordField: string = 'password'): Rule => ({
    validator: (_, value) => {
      if (!value || getFieldValue(passwordField) === value) {
        return Promise.resolve()
      }
      return Promise.reject(new Error('The two passwords do not match!'))
    },
  }),

  /**
   * Phone number validation (flexible format)
   * Accepts various formats: (123) 456-7890, 123-456-7890, 1234567890, etc.
   */
  phone: {
    pattern: /^[\d\s()+-]+$/,
    message: 'Please input a valid phone number!',
  },

  /**
   * Phone number validation (strict US format)
   * Format: (123) 456-7890 or 123-456-7890
   */
  phoneStrict: {
    pattern: /^(\(\d{3}\)\s?|\d{3}-)\d{3}-\d{4}$/,
    message: 'Please input a valid phone number (e.g., (123) 456-7890)!',
  },

  /**
   * Alphanumeric validation (letters and numbers only)
   */
  alphanumeric: {
    pattern: /^[a-zA-Z0-9]+$/,
    message: 'Only letters and numbers are allowed!',
  },

  /**
   * Alpha validation (letters only)
   */
  alpha: {
    pattern: /^[a-zA-Z]+$/,
    message: 'Only letters are allowed!',
  },

  /**
   * Numeric validation (numbers only)
   */
  numeric: {
    pattern: /^[0-9]+$/,
    message: 'Only numbers are allowed!',
  },

  /**
   * Decimal number validation
   *
   * @param decimals - Maximum decimal places (optional)
   * @returns Validation rule
   */
  decimal: (decimals?: number): Rule => {
    const pattern = decimals ? new RegExp(`^\\d+(\\.\\d{1,${decimals}})?$`) : /^\d+(\.\d+)?$/

    return {
      pattern,
      message: decimals
        ? `Must be a number with at most ${decimals} decimal places!`
        : 'Must be a valid decimal number!',
    }
  },

  /**
   * Positive number validation
   */
  positiveNumber: {
    type: 'number' as const,
    min: 0,
    message: 'Must be a positive number!',
  },

  /**
   * Integer validation
   */
  integer: {
    type: 'integer' as const,
    message: 'Must be an integer!',
  },

  /**
   * Array minimum length validation
   *
   * @param min - Minimum number of items
   * @returns Validation rule
   */
  arrayMinLength: (min: number): Rule => ({
    type: 'array' as const,
    min,
    message: `Please select at least ${min} item${min > 1 ? 's' : ''}!`,
  }),

  /**
   * Array maximum length validation
   *
   * @param max - Maximum number of items
   * @returns Validation rule
   */
  arrayMaxLength: (max: number): Rule => ({
    type: 'array' as const,
    max,
    message: `Please select at most ${max} item${max > 1 ? 's' : ''}!`,
  }),

  /**
   * Array required validation
   */
  arrayRequired: {
    type: 'array' as const,
    required: true,
    message: 'Please select at least one item!',
  },

  /**
   * Date required validation
   */
  dateRequired: {
    type: 'date' as const,
    required: true,
    message: 'Please select a date!',
  },

  /**
   * Future date validation
   */
  futureDate: {
    type: 'date' as const,
    validator: (_: unknown, value: unknown) => {
      if (!value || value > new Date()) {
        return Promise.resolve()
      }
      return Promise.reject(new Error('Date must be in the future!'))
    },
  },

  /**
   * Past date validation
   */
  pastDate: {
    type: 'date' as const,
    validator: (_: unknown, value: unknown) => {
      if (!value || value < new Date()) {
        return Promise.resolve()
      }
      return Promise.reject(new Error('Date must be in the past!'))
    },
  },

  /**
   * JSON validation
   */
  json: {
    validator: (_: unknown, value: unknown) => {
      if (!value) {
        return Promise.resolve()
      }
      try {
        JSON.parse(value as string)
        return Promise.resolve()
      } catch {
        return Promise.reject(new Error('Must be valid JSON!'))
      }
    },
  },

  /**
   * Whitespace not allowed validation
   */
  noWhitespace: {
    pattern: /^\S+$/,
    message: 'Whitespace is not allowed!',
  },

  /**
   * Slug validation (URL-friendly string)
   * Lowercase letters, numbers, and hyphens only
   */
  slug: {
    pattern: /^[a-z0-9-]+$/,
    message: 'Only lowercase letters, numbers, and hyphens are allowed!',
  },

  /**
   * HEX color validation
   */
  hexColor: {
    pattern: /^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$/,
    message: 'Please input a valid HEX color (e.g., #FF5733)!',
  },

  /**
   * IPv4 address validation
   */
  ipv4: {
    pattern: /^((25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.){3}(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)$/,
    message: 'Please input a valid IPv4 address!',
  },

  /**
   * Coordinate validation (latitude)
   */
  latitude: {
    type: 'number' as const,
    min: -90,
    max: 90,
    message: 'Latitude must be between -90 and 90!',
  },

  /**
   * Coordinate validation (longitude)
   */
  longitude: {
    type: 'number' as const,
    min: -180,
    max: 180,
    message: 'Longitude must be between -180 and 180!',
  },

  /**
   * File size validation
   *
   * @param maxSizeMB - Maximum file size in MB
   * @returns Validation rule
   *
   * @example
   * validationRules.fileSize(10) // Max 10MB
   */
  fileSize: (maxSizeMB: number): Rule => ({
    validator: (_, file) => {
      if (!file) {
        return Promise.resolve()
      }

      const fileSizeMB = file.size / 1024 / 1024

      if (fileSizeMB <= maxSizeMB) {
        return Promise.resolve()
      }

      return Promise.reject(new Error(`File size must be less than ${maxSizeMB}MB!`))
    },
  }),

  /**
   * File type validation
   *
   * @param allowedTypes - Array of allowed MIME types or extensions
   * @returns Validation rule
   *
   * @example
   * validationRules.fileType(['image/jpeg', 'image/png', '.jpg', '.png'])
   */
  fileType: (allowedTypes: string[]): Rule => ({
    validator: (_, file) => {
      if (!file) {
        return Promise.resolve()
      }

      const fileType = file.type
      const fileName = file.name
      const fileExtension = fileName.substring(fileName.lastIndexOf('.')).toLowerCase()

      const isValid =
        allowedTypes.some((type) => type === fileType) ||
        allowedTypes.some((type) => type.startsWith('.') && type === fileExtension)

      if (isValid) {
        return Promise.resolve()
      }

      return Promise.reject(new Error(`Only ${allowedTypes.join(', ')} files are allowed!`))
    },
  }),

  /**
   * Custom async validator
   * Useful for server-side validation (e.g., check if username exists)
   *
   * @param validateFn - Async validation function that returns boolean or throws error
   * @param errorMessage - Error message if validation fails
   * @returns Validation rule
   *
   * @example
   * validationRules.asyncValidator(
   *   async (value) => {
   *     const exists = await checkUsernameExists(value)
   *     return !exists
   *   },
   *   'Username already exists!'
   * )
   */
  asyncValidator: (validateFn: (value: any) => Promise<boolean>, errorMessage: string): Rule => ({
    validator: async (_: unknown, value) => {
      if (!value) {
        return Promise.resolve()
      }

      try {
        const isValid = await validateFn(value)
        if (isValid) {
          return Promise.resolve()
        }
        return Promise.reject(new Error(errorMessage))
      } catch (error) {
        return Promise.reject(new Error(errorMessage))
      }
    },
  }),
}

/**
 * Common field validation presets
 * Combines multiple rules for common field types
 */
export const fieldPresets = {
  /**
   * Full name validation
   * - Required
   * - 2-100 characters
   * - Letters, spaces, hyphens, and apostrophes only
   */
  fullName: [
    validationRules.required('full name'),
    validationRules.lengthRange(2, 100),
    validationRules.pattern(/^[a-zA-Z\s'-]+$/, 'Only letters, spaces, hyphens, and apostrophes are allowed!'),
  ],

  /**
   * Organization name validation
   * - 2-200 characters
   */
  organization: [validationRules.lengthRange(2, 200)],

  /**
   * Project name validation
   * - Required
   * - 3-100 characters
   */
  projectName: [validationRules.required('project name'), validationRules.lengthRange(3, 100)],

  /**
   * Sample name validation
   * - Required
   * - 1-100 characters
   */
  sampleName: [validationRules.required('sample name'), validationRules.lengthRange(1, 100)],

  /**
   * Description validation
   * - Optional
   * - Maximum 1000 characters
   */
  description: [validationRules.maxLength(1000)],

  /**
   * Search query validation
   * - Maximum 200 characters
   */
  searchQuery: [validationRules.maxLength(200)],
}

export default validationRules
